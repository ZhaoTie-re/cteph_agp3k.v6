#!/usr/bin/env python3
"""Step 04 — one confidence tier per reported locus.

Three criteria (config.py records which thresholds were fixed before step 03 and
which were revised after it, and why):

  unconfirmed   share of control chromosomes on an allele ToMMo does not carry although
                its dictionary could (named by IPD-IMGT no later than ToMMo's newest
                allele) <= TIER_MAX_UNCONFIRMED (ToMMo loci)
  dropout       lower 95 % bound of the share of heterozygotes reported as homozygotes
                <= TIER_MAX_DROPOUT (obligate genes; Hardy-Weinberg, own frequencies)
  call quality  share of high-confidence calls >= TIER_MIN_HIGH_CONF (every locus)

Frequency discordance against ToMMo and 1KG JPT is carried along, not tiered.

Tier = 1 + number of criteria failed, capped at 3. A locus with no external panel
starts at 2, because nothing outside our own data speaks for it.

Output: results/04.confidence/locus_tiers.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def main():
    loc = pd.read_csv(C.D_PANELS / 'locus_concordance.tsv', sep='\t')
    hz = pd.read_csv(C.D_PANELS / 'homozygosity.tsv', sep='\t').set_index('gene')
    gs = pd.read_csv(C.D_EVIDENCE / 'gene_summary.tsv', sep='\t').set_index('gene')
    tom = loc[loc.panel == 'ToMMo'].set_index('gene')
    kg = loc[loc.panel == '1KG JPT'].set_index('gene')

    rows = []
    for g in C.REPORTED_GENES:
        r = {'gene': g, 'external_panel': 'ToMMo' if g in tom.index else 'none'}
        crit = {}
        if g in tom.index and tom.at[g, 'frequency_comparable'] == 1:
            r['excess_D'] = tom.at[g, 'excess_D']
            r['D_shared'] = tom.at[g, 'D_shared']
        if g in tom.index:
            r['unconfirmed_any'] = tom.at[g, 'unconfirmed']
            r['unconfirmed_after_vintage'] = tom.at[g, 'unconfirmed_after_vintage']
            r['unconfirmed'] = tom.at[g, 'unconfirmed_assessable']
            r['unconfirmed_high'] = tom.at[g, 'unconfirmed_assessable_high']
            crit['unconfirmed'] = r['unconfirmed'] <= C.TIER_MAX_UNCONFIRMED
        if g in hz.index:
            r['hom_ratio'] = hz.at[g, 'ratio']
            r['dropout'], r['dropout_lo'], r['dropout_hi'] = (hz.at[g, 'dropout'], hz.at[g, 'dropout_lo'],
                                                              hz.at[g, 'dropout_hi'])
            crit['dropout'] = r['dropout_lo'] <= C.TIER_MAX_DROPOUT
        r['share_high'] = gs.at[g, 'share_high']
        crit['call_quality'] = r['share_high'] >= C.TIER_MIN_HIGH_CONF
        if g in kg.index:
            r['excess_D_1kg'] = kg.at[g, 'excess_D']
        for k in ('unconfirmed', 'dropout', 'call_quality'):
            r[f'pass_{k}'] = '' if k not in crit else int(crit[k])
        failed = [k for k, v in crit.items() if not v]
        start = 1 if r['external_panel'] != 'none' else 2
        r['criteria_applied'] = len(crit)
        r['criteria_failed'] = ';'.join(failed)
        r['tier'] = min(start + len(failed), 3)
        rows.append(r)
    out = pd.DataFrame(rows)
    C.D_CONF.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.D_CONF / 'locus_tiers.tsv', sep='\t', index=False, float_format='%.6g')
    pd.set_option('display.width', 220)
    print(out.round(4).to_string(index=False))
    print('[tiers]', out.tier.value_counts().sort_index().to_dict())


if __name__ == '__main__':
    main()
