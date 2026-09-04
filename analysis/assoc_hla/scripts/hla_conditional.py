#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Forward stepwise conditional analysis over the MHC -- how many
#           INDEPENDENT signals there are, not how many markers are significant.
#
#           HLA markers are not independent tests and a single causal allele
#           drags a whole haplotype over the threshold with it. The question a
#           reader actually has is "how many separate things are going on here",
#           and the accepted way to answer it is to add the lead marker's dosage
#           to the covariate set and re-scan, repeating until nothing clears the
#           threshold (Hirata et al. Nat. Genet. 2019; Nat. Protoc. 2023).
#
#           CONDITIONING IS GLOBAL ACROSS THE MHC, not per gene. LD in this
#           region spans genes -- DRB1, DQA1 and DQB1 travel together on
#           haplotypes -- so conditioning within a gene would leave the signal
#           standing in its neighbour and count one thing twice.
#
#           THE NULL MODEL IS NEVER REFITTED under the random model. Conditioning
#           adds a covariate at the association step only, exactly as in
#           assoc_saige, so the rounds stay comparable with one another and with
#           round 0.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import os
import subprocess
import sys

import pandas as pd

CANON_MIN = ['ID', '#CHROM', 'POS', 'A1', 'A1_FREQ', 'OR', 'LOG(OR)_SE', 'P', 'ERRCODE']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bfile', required=True)
    p.add_argument('--extract', required=True, help='the tested marker list for this class')
    p.add_argument('--marker-map', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True, choices=('fixed', 'random'))
    p.add_argument('--marker-class', required=True)
    p.add_argument('--threshold', type=float, required=True,
                   help='a round passes while its top P is below this')
    p.add_argument('--max-rounds', type=int, default=5)
    p.add_argument('--max-frac-unconfirmed', type=float, default=0.5,
                   help='a residue whose MINOR-state chromosomes come mostly from alleles '
                        'the reference panel has never seen is not eligible to LEAD a '
                        'round. It is still tested and still appears in every round table; '
                        'it just cannot be the thing the next round conditions on.')
    p.add_argument('--max-frac-unconfirmed-allele-flag', dest='drop_absent_alleles',
                   action='store_true', default=True,
                   help='same rule for allele markers, via in_reference == 0')
    p.add_argument('--adapter', required=True, help='hla_to_sumstats.py')
    # fixed
    p.add_argument('--plink2'), p.add_argument('--pheno'), p.add_argument('--covar')
    p.add_argument('--covar-name'), p.add_argument('--pheno-name', default='PHENO1')
    # random
    p.add_argument('--saige-step2'), p.add_argument('--saige-rscript')
    p.add_argument('--gmmat'), p.add_argument('--variance-ratio')
    p.add_argument('--p-value', default='spa')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out-prefix', required=True)
    return p.parse_args()


def run(cmd, log):
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        tail = open(log).read()[-1500:]
        raise SystemExit(f'ABORT: command failed ({r.returncode})\n  {" ".join(cmd)}\n{tail}')


def scan(args, rnd, conditioned):
    """One association pass, conditioned on `conditioned`. Returns the canonical table."""
    tag = f'round{rnd}'
    if args.model == 'fixed':
        cmd = [args.plink2, '--bfile', args.bfile, '--extract', args.extract,
               '--pheno', args.pheno, '--pheno-name', args.pheno_name, '--1',
               '--covar', args.covar, '--covar-name', args.covar_name,
               '--glm', 'omit-ref', 'hide-covar', 'firth-fallback',
               '--ci', '0.95', '--threads', str(args.threads), '--out', tag]
        if conditioned:
            with open(f'{tag}.cond', 'w') as fh:
                fh.write('\n'.join(conditioned) + '\n')
            cmd += ['--condition-list', f'{tag}.cond']
        run(cmd, f'{tag}.plink.log')
        hit = [f for f in os.listdir('.')
               if f.startswith(tag) and f.endswith('.glm.logistic.hybrid')]
        if not hit:
            raise SystemExit(f'ABORT: plink2 wrote no .glm.logistic.hybrid for {tag}')
        raw, engine = hit[0], 'plink2'
    else:
        cmd = [args.saige_rscript, args.saige_step2,
               f'--bedFile={args.bfile}.bed', f'--bimFile={args.bfile}.bim',
               f'--famFile={args.bfile}.fam', '--AlleleOrder=alt-first', '--chrom=6',
               f'--GMMATmodelFile={args.gmmat}',
               f'--varianceRatioFile={args.variance_ratio}',
               f'--SAIGEOutputFile={tag}.assoc.txt', '--is_Firth_beta=TRUE',
               '--LOCO=TRUE', '--is_output_moreDetails=TRUE']
        if conditioned:
            cmd.append('--condition=' + ','.join(conditioned))
        run(cmd, f'{tag}.saige.log')
        raw, engine = f'{tag}.assoc.txt', 'saige'

    out = f'{tag}.sumstat.tsv'
    run([sys.executable, args.adapter, '--assoc', raw, '--engine', engine,
         '--bim', f'{args.bfile}.bim', '--marker-map', args.marker_map,
         '--cohort', args.cohort, '--model', args.model,
         '--marker-class', args.marker_class, '--p-value', args.p_value,
         '--out-sumstats', out, '--out-extra', f'{tag}.extra.tsv'], f'{tag}.adapt.log')
    d = pd.read_csv(out, sep='\t')
    return d[d.ERRCODE == '.'].copy(), out


def main():
    args = parse_args()
    locus = f'mhc_{args.marker_class}'

    # CONDITIONING ON AN ARTEFACT REORDERS EVERYTHING BELOW IT.
    # The cascade decides which signals are INDEPENDENT, so whatever leads round 0
    # is conditioned out of every later round. If that lead is a typing artefact,
    # the real signals are pushed down the list and reported as secondary.
    # Measured on this cohort before the guard: round 0 was led by AA_A_178_T at
    # P = 7e-9, whose frac_unconfirmed_allele is 0.85 -- 15 of its 20 minor-state
    # chromosomes carry one allele a 61,424-person Japanese panel has never seen.
    # The pipeline had already computed that flag and this script ignored it.
    #
    # Flagged markers are NOT removed: they are tested every round and appear in
    # every round table. They are only barred from being the lead.
    qc = pd.read_csv(args.marker_map, sep='\t')
    ineligible = set()
    if 'frac_unconfirmed_allele' in qc.columns:
        bad = pd.to_numeric(qc.frac_unconfirmed_allele, errors='coerce')
        ineligible |= set(qc.loc[bad > args.max_frac_unconfirmed, 'id'])
    if args.drop_absent_alleles and 'in_reference' in qc.columns:
        ref = pd.to_numeric(qc.in_reference, errors='coerce')
        ineligible |= set(qc.loc[ref == 0, 'id'])

    conditioned, rounds, signals = [], [], []

    for rnd in range(args.max_rounds + 1):
        d, path = scan(args, rnd, conditioned)
        d = d[~d.ID.isin(conditioned)]
        if not len(d):
            rounds.append({'round': rnd, 'n_conditioned_on': len(conditioned),
                           'conditioned_on': ';'.join(conditioned) or 'none',
                           'n_variants': 0, 'top_id': 'NA', 'top_p': float('nan'),
                           'top_or': float('nan'), 'passes_threshold': 0,
                           'cohort': args.cohort, 'locus_id': locus})
            break
        pool = d[~d.ID.isin(ineligible)]
        if not len(pool):
            rounds.append({'round': rnd, 'n_conditioned_on': len(conditioned),
                           'conditioned_on': ';'.join(conditioned) or 'none',
                           'n_variants': len(d), 'top_id': 'NA', 'top_p': float('nan'),
                           'top_or': float('nan'), 'passes_threshold': 0,
                           'cohort': args.cohort, 'locus_id': locus})
            break
        top = pool.loc[pool.P.idxmin()]
        passes = bool(top.P < args.threshold)
        rounds.append({'round': rnd, 'n_conditioned_on': len(conditioned),
                       'conditioned_on': ';'.join(conditioned) or 'none',
                       'n_variants': len(d), 'top_id': top.ID, 'top_p': float(top.P),
                       'top_or': float(top.OR), 'passes_threshold': int(passes),
                       'cohort': args.cohort, 'locus_id': locus})
        os.replace(path, f'{args.out_prefix}.round{rnd}.tsv')
        if not passes:
            break
        signals.append({'cohort': args.cohort, 'locus_id': locus, 'signal': len(signals) + 1,
                        'variant_id': top.ID, 'chrom': top['#CHROM'], 'pos': top.POS,
                        'a1': top.A1, 'a1_freq': top.A1_FREQ, 'p': float(top.P),
                        'or': float(top.OR), 'se': float(top['LOG(OR)_SE']),
                        'conditioned_on': ';'.join(conditioned) or 'none'})
        conditioned.append(top.ID)

    pd.DataFrame(rounds).to_csv(f'{args.out_prefix}.rounds.tsv', sep='\t',
                                index=False, na_rep='NA')
    pd.DataFrame(signals, columns=['cohort', 'locus_id', 'signal', 'variant_id', 'chrom',
                                   'pos', 'a1', 'a1_freq', 'p', 'or', 'se', 'conditioned_on']
                 ).to_csv(f'{args.out_prefix}.signals.tsv', sep='\t', index=False, na_rep='NA')

    print(f'[hla_conditional] {args.cohort}/{args.model}/{args.marker_class}: '
          f'{len(signals)} independent signal(s) at P < {args.threshold:.3e} '
          f'over {len(rounds)} round(s)')
    if ineligible:
        print(f'    {len(ineligible)} marker(s) barred from leading a round on the '
              f'artefact flags; they are still tested and still in every round table')
    for s in signals:
        print(f'    signal {s["signal"]}: {s["variant_id"]}  P={s["p"]:.3e}  OR={s["or"]:.3f}')
    if len(signals) == args.max_rounds:
        print(f'    NOTE: stopped at --max-rounds {args.max_rounds}, not at the threshold; '
              f'there may be more signals.', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in hla_conditional: {e}', file=sys.stderr)
        sys.exit(1)
