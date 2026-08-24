#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Stepwise conditional analysis inside one genome-wide locus, under the
#           mixed model. Round 0 is the unconditioned fit restricted to the locus
#           window; each later round adds the previous round's top variant to
#           SAIGE's --condition and re-fits. The loop stops when nothing in the
#           window still passes the genome-wide threshold, and the number of
#           rounds that yielded a signal is the number of independent signals.
#
#           Same contract as the fixed-model component's conditional analysis —
#           the same three output files, the same columns — so the shared
#           conditional figure reads either engine's output unchanged.
#
#           TWO THINGS ARE SPECIFIC TO THE MIXED MODEL:
#
#           1. THE NULL MODEL IS NEVER RE-FITTED. Conditioning changes which
#              variants are covariates in step 2 only; the GRM and the variance
#              ratio from step 1 are re-used across every round. Re-fitting would
#              cost hours per round and would also change the null the rounds are
#              compared against, making them incomparable.
#
#           2. THE CONDITIONING ID IS SAIGE'S OWN. SAIGE names variants by its
#              own convention, which need not match the .bim or our chr:pos:REF:ALT
#              key, so the id fed back to --condition is taken verbatim from the
#              MarkerID SAIGE printed for that variant. Self-consistent by
#              construction, rather than by assuming a format.
#
#           The window is cut out with plink2 first so each round tests a few
#           thousand variants rather than a whole chromosome.
# Component: assoc_saige
# Used by : assoc_saige.nf  process CONDITIONAL
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from saige_to_sumstats import read_saige, to_canonical


def parse_args():
    p = argparse.ArgumentParser(description='Stepwise conditional analysis for one locus (SAIGE).')
    p.add_argument('--saige-rscript', required=True)
    p.add_argument('--saige-step2', required=True)
    p.add_argument('--plink2', default='plink2')
    p.add_argument('--bed', required=True)
    p.add_argument('--bim', required=True)
    p.add_argument('--fam', required=True)
    p.add_argument('--rda', required=True, help='step-1 null model (.rda)')
    p.add_argument('--variance-ratio', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--locus-id', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--end', type=int, required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--p-value', default='spa', choices=['spa', 'no_spa'])
    p.add_argument('--p-threshold', type=float, default=5e-8)
    p.add_argument('--max-rounds', type=int, default=5)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out-dir', required=True)
    return p.parse_args()


def cut_window(args, work):
    """The locus window as a small bfile, so a round is seconds rather than a
    whole-chromosome pass."""
    pfx = work / 'window'
    cmd = [args.plink2, '--bed', args.bed, '--bim', args.bim, '--fam', args.fam,
           '--chr', str(args.chrom).replace('chr', ''),
           '--from-bp', str(args.start), '--to-bp', str(args.end),
           '--make-bed', '--out', str(pfx), '--threads', str(args.threads)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not pfx.with_suffix('.bed').exists():
        raise SystemExit(f'plink2 produced no window for {args.locus_id}\n{r.stderr[-2000:]}')
    return pfx


def run_step2(args, window, out_file, condition_ids):
    """One SAIGE step-2 pass over the window, optionally conditioned."""
    cmd = [args.saige_rscript, args.saige_step2,
           f'--bedFile={window}.bed', f'--bimFile={window}.bim', f'--famFile={window}.fam',
           '--AlleleOrder=alt-first',
           f'--chrom={str(args.chrom).replace("chr", "")}',
           f'--GMMATmodelFile={args.rda}',
           f'--varianceRatioFile={args.variance_ratio}',
           f'--SAIGEOutputFile={out_file}',
           '--is_Firth_beta=FALSE', '--LOCO=TRUE', '--is_output_moreDetails=TRUE']
    if condition_ids:
        cmd.append('--condition=' + ','.join(condition_ids))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if not Path(out_file).exists():
        # Not a crash: with every variant in the window collinear with the
        # conditioning set SAIGE has nothing left to test. That is the outcome
        # the stopping rule is looking for.
        print(f'[saige_conditional] no output for {out_file}\n{res.stderr[-2000:]}', file=sys.stderr)
        return None
    d = read_saige(out_file)
    if not len(d):
        return None
    # WHEN CONDITIONING, THE ANSWER IS IN THE `_c` COLUMNS. SAIGE reports the
    # MARGINAL statistics in BETA/SE/p.value even when --condition is given, and
    # puts the conditional ones in BETA_c/SE_c/p.value_c. Reading the plain
    # columns therefore returns the unconditioned result every round: the loop
    # never converges, re-selects the same variant, and reports one "independent
    # signal" per round. Verified against SAIGE's own output: conditioning on a
    # variant drives that variant's p.value_c to ~1 while p.value is unchanged.
    if condition_ids:
        pairs = [('BETA_c', 'BETA'), ('SE_c', 'SE'), ('Tstat_c', 'Tstat'),
                 ('var_c', 'var'), ('p.value_c', 'p.value'),
                 ('p.value.NA_c', 'p.value.NA')]
        missing = [c for c, _ in pairs if c not in d.columns]
        if missing:
            raise SystemExit(f'conditioned round produced no {missing} columns — SAIGE did not '
                             f'apply --condition; check the marker id format')
        for src, dst in pairs:
            d[dst] = d[src]
    canon, _pcol, _bad = to_canonical(d, args.p_value)
    canon = canon.rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
    canon['MarkerID_saige'] = d['MarkerID'].astype(str).values if 'MarkerID' in d else canon['ID']
    return canon[canon['ERRCODE'].eq('.') & canon['P'].notna()]


def main():
    args = parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    work = out / '_work'; work.mkdir(exist_ok=True)
    window = cut_window(args, work)

    signals, rounds, condition = [], [], []
    for rnd in range(args.max_rounds + 1):
        df = run_step2(args, window, work / f'{args.locus_id}.round{rnd}.saige.txt', condition)
        if df is None or not len(df):
            rounds.append({'round': rnd, 'n_conditioned_on': len(condition),
                           'conditioned_on': ','.join(condition) or '.',
                           'n_variants': 0, 'top_id': '.', 'top_p': None,
                           'top_or': None, 'passes_threshold': False})
            break
        top = df.loc[df['P'].idxmin()]
        # A conditioned round that picks a variant already in the conditioning
        # set means --condition did nothing. That must fail loudly: silently it
        # looks like one independent signal per round, forever.
        if condition and str(top['MarkerID_saige']) in condition:
            raise SystemExit(
                f"round {rnd} re-selected {top['MarkerID_saige']}, which is already conditioned "
                f"on — the conditioning had no effect and the result would be meaningless")
        passes = bool(top['P'] < args.p_threshold)
        rounds.append({'round': rnd, 'n_conditioned_on': len(condition),
                       'conditioned_on': ','.join(condition) or '.',
                       'n_variants': len(df), 'top_id': top['ID'],
                       'top_p': top['P'], 'top_or': top.get('OR'),
                       'passes_threshold': passes})
        df.drop(columns=['MarkerID_saige']) \
          .assign(round=rnd, conditioned_on=','.join(condition) or '.') \
          .to_csv(out / f'{args.locus_id}.round{rnd}.tsv', sep='\t', index=False)
        if not passes:
            break
        signals.append({'cohort': args.cohort, 'locus_id': args.locus_id, 'signal': len(signals) + 1,
                        'variant_id': top['ID'], 'chrom': top['CHROM'], 'pos': int(top['POS']),
                        'a1': top.get('A1', ''), 'a1_freq': top.get('A1_FREQ'),
                        'p': top['P'], 'or': top.get('OR'), 'se': top.get('SE'),
                        'conditioned_on': ','.join(condition) or '.'})
        # SAIGE's own name for the variant — see the header note.
        condition.append(str(top['MarkerID_saige']))

    pd.DataFrame(rounds).assign(cohort=args.cohort, locus_id=args.locus_id) \
      .to_csv(out / f'{args.locus_id}.rounds.tsv', sep='\t', index=False)
    pd.DataFrame(signals, columns=['cohort', 'locus_id', 'signal', 'variant_id', 'chrom', 'pos',
                                   'a1', 'a1_freq', 'p', 'or', 'se', 'conditioned_on']) \
      .to_csv(out / f'{args.locus_id}.signals.tsv', sep='\t', index=False)

    print(f'[saige_conditional] {args.cohort} {args.locus_id}: {len(signals)} independent '
          f'signal(s) over {len(rounds)} round(s); lead {args.lead_id}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in saige_conditional: {e}', file=sys.stderr)
        sys.exit(1)
