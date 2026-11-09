#!/usr/bin/env python3
"""Step 02 — the read-level evidence behind each call, and a confidence per call.

Thresholds come from config.py, where they were fixed before this step was run.

Per call (reported genes):
  low_support     an allele with fewer than MIN_ALLELE_READS reads assigned to it
  imbalanced      a heterozygote whose minor allele has less than MIN_MINOR_SHARE of the
                  two alleles' UNIQUE reads (reads assigned to one allele and not the
                  other). Assessed only with at least MIN_UNIQUE_READS unique reads:
                  closely related alleles differ at few positions and leave too few.
  drb_*           from step 01's DRB1-based copy-number reconciliation
  2-field basis   agreed (reads) > prior with ratio >= PRIOR_DOMINANCE > P group /
                  weaker prior

Confidence:
  high     2-field alleles agreed by the reads, or chosen by HLA-HD's prior over
           alternatives at least PRIOR_DOMINANCE times rarer; and no flag
  medium   one weakness: a P-group or weak-prior resolution, or one flag
  low      two or more weaknesses, or no 2-field allele (ambiguous)
  no_call  failed, unresolved, or a DRB3/4/5 copy expected but not typed

Per sample: read pairs, failed obligate genes, reported genes without a 2-field
allele, low-confidence calls. Outliers (listed, not removed): read pairs below the
OUTLIER_READ_PAIRS_PCTL percentile, any failed obligate gene, or at least
OUTLIER_MIN_UNRESOLVED reported genes without a 2-field allele.

Outputs (results/02.evidence/): call_evidence.tsv, gene_summary.tsv, sample_qc.tsv,
drb_consistency.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def weaknesses(r):
    """Named weaknesses of one call; None when there is no call."""
    if r.state in C.NO_CALL_STATES:
        return None
    if r.state == 'absent':
        return [r.drb_flag] if r.drb_flag else []
    w = ['no_2field'] if r.state == 'ambiguous' else []
    for k in (('1',) if r.copies == 1 else ('1', '2')):
        b = getattr(r, f'basis_{k}')
        if b == 'pgroup':
            w.append('pgroup')
        elif b == 'prior':
            q = pd.to_numeric(getattr(r, f'prior_ratio_{k}'), errors='coerce')
            if not (q >= C.PRIOR_DOMINANCE):
                w.append('weak_prior')
    if r.low_support:
        w.append('low_support')
    if r.imbalanced:
        w.append('imbalanced')
    if r.drb_flag:
        w.append(r.drb_flag)
    return list(dict.fromkeys(w))


def confidence(w):
    if w is None:
        return 'no_call'
    if 'no_2field' in w:
        return 'low'
    return 'high' if not w else 'medium' if len(w) == 1 else 'low'


def main():
    calls = pd.read_csv(C.D_CALLS / 'calls_long.tsv', sep='\t', dtype=str).fillna('')
    samples = pd.read_csv(C.D_INPUTS / 'sample_inputs.tsv', sep='\t', dtype=str)
    rep = calls[calls.reported == '1'].copy()
    for c in ('reads_1', 'reads_2', 'unique_1', 'unique_2', 'shared', 'copies'):
        rep[c] = pd.to_numeric(rep[c], errors='coerce')

    called = ~rep.state.isin(C.NO_CALL_STATES + ['absent'])
    # Paralogue calls reduced to one copy: support is that copy's reads only.
    one_copy = rep.copies == 1
    reads_min = np.where(one_copy, rep.reads_1, rep[['reads_1', 'reads_2']].min(axis=1))
    rep['min_allele_reads'] = reads_min
    rep['low_support'] = called & (pd.Series(reads_min, index=rep.index).fillna(0) < C.MIN_ALLELE_READS)

    het = (rep.state == 'heterozygous') | ((rep.state == 'present') & (rep.copies == 2)
                                           & (rep.allele_1 != rep.allele_2))
    uniq = rep.unique_1.fillna(0) + rep.unique_2.fillna(0)
    rep['unique_total'] = np.where(het, uniq, np.nan)
    rep['minor_share'] = np.where(het & (uniq > 0),
                                  rep[['unique_1', 'unique_2']].min(axis=1) / uniq.replace(0, np.nan),
                                  np.nan)
    rep['balance_assessed'] = het & (uniq >= C.MIN_UNIQUE_READS)
    rep['imbalanced'] = rep.balance_assessed & (rep.minor_share < C.MIN_MINOR_SHARE)

    ws = [weaknesses(r) for r in rep.itertuples(index=False)]
    rep['weaknesses'] = [';'.join(w) if w else '' for w in ws]
    rep['confidence'] = [confidence(w) for w in ws]

    C.D_EVIDENCE.mkdir(parents=True, exist_ok=True)
    cols = ['sample_id', 'gene', 'state', 'copies', 'allele_1', 'allele_2', 'basis_1', 'basis_2',
            'prior_ratio_1', 'prior_ratio_2', 'reads_1', 'reads_2', 'unique_1', 'unique_2',
            'shared', 'min_allele_reads', 'unique_total', 'minor_share', 'balance_assessed',
            'low_support', 'imbalanced', 'drb_flag', 'weaknesses', 'confidence']
    rep[cols].to_csv(C.D_EVIDENCE / 'call_evidence.tsv', sep='\t', index=False)

    # ---- per gene
    g = []
    for gene in C.REPORTED_GENES:
        s = rep[rep.gene == gene]
        n = len(s)
        alleles = pd.concat([s[['basis_1', 'allele_1']].set_axis(['b', 'a'], axis=1),
                             s[['basis_2', 'allele_2']].set_axis(['b', 'a'], axis=1)])
        alleles = alleles[alleles.a != '']
        hs = s[s.balance_assessed]
        g.append({
            'gene': gene, 'class': ('paralogue' if gene in C.PARALOGUES else
                                    'class I' if gene in C.CLASS_I else 'class II'),
            'samples': n,
            **{f'state_{k}': int((s.state == k).sum()) for k in
               ('heterozygous', 'homozygous', 'present', 'absent', 'ambiguous',
                'unresolved', 'failed', 'missing')},
            'alleles_2field': len(alleles),
            **{f'basis_{k}': int((alleles.b == k).sum()) for k in ('agreed', 'prior', 'pgroup')},
            'weak_prior': int(rep[rep.gene == gene].weaknesses.str.contains('weak_prior').sum()),
            'median_min_reads': float(np.nanmedian(s.min_allele_reads)) if s.min_allele_reads.notna().any() else np.nan,
            'low_support': int(s.low_support.sum()),
            'het_assessed': len(hs), 'imbalanced': int(s.imbalanced.sum()),
            **{f'conf_{k}': int((s.confidence == k).sum()) for k in ('high', 'medium', 'low', 'no_call')},
        })
    gs = pd.DataFrame(g)
    gs['share_high'] = (gs.conf_high / gs.samples).round(4)
    gs.to_csv(C.D_EVIDENCE / 'gene_summary.tsv', sep='\t', index=False)

    # ---- DRB1 vs HLA-HD paralogue copy number, before reconciliation
    rec = pd.read_csv(C.D_CALLS / 'drb_reconciliation.tsv', sep='\t', dtype=str).fillna('')
    rec['expected'] = rec.expected.replace('', 'unknown')
    tab = rec.groupby(['gene', 'expected', 'observed']).size().rename('samples').reset_index()
    tab.to_csv(C.D_EVIDENCE / 'drb_consistency.tsv', sep='\t', index=False)

    # ---- per sample
    rep['no2'] = rep.state.isin(['ambiguous', 'unresolved'])
    per = rep.groupby('sample_id').agg(
        failed_obligate=('state', lambda v: int((v == 'failed').sum())),
        no_2field=('no2', 'sum'),
        low_conf=('confidence', lambda v: int(v.isin(['low', 'no_call']).sum())),
        high_conf=('confidence', lambda v: int((v == 'high').sum())),
        flags_imbalanced=('imbalanced', 'sum'), flags_low_support=('low_support', 'sum'))
    sq = samples.set_index('sample_id')[['group', 'sex', 'platform', 'observed_depth',
                                         'read_pairs']].join(per)
    sq['read_pairs'] = sq.read_pairs.astype(int)
    cut = np.percentile(sq.read_pairs, C.OUTLIER_READ_PAIRS_PCTL)
    sq['low_reads'] = sq.read_pairs < cut
    sq['outlier'] = (sq.low_reads | (sq.failed_obligate > 0)
                     | (sq.no_2field >= C.OUTLIER_MIN_UNRESOLVED))
    sq['outlier_reason'] = [';'.join(x for x, f in (('read_pairs<p1', a), ('failed_obligate', b > 0),
                                                       ('no_2field>=3', c >= C.OUTLIER_MIN_UNRESOLVED)) if f)
                            for a, b, c in zip(sq.low_reads, sq.failed_obligate, sq.no_2field)]
    sq.reset_index().to_csv(C.D_EVIDENCE / 'sample_qc.tsv', sep='\t', index=False)

    print(gs[['gene', 'samples', 'state_heterozygous', 'state_homozygous', 'state_present',
              'state_absent', 'state_ambiguous', 'state_failed', 'state_missing', 'basis_prior',
              'low_support', 'het_assessed', 'imbalanced', 'conf_high', 'conf_medium',
              'conf_low', 'conf_no_call', 'share_high']].to_string(index=False))
    print(f'[evidence] read-pair cut (p{C.OUTLIER_READ_PAIRS_PCTL:g}) = {cut:,.0f}; '
          f'outlier samples: {int(sq.outlier.sum())} '
          f'({sq.outlier_reason[sq.outlier].str.split(";").explode().value_counts().to_dict()})')


if __name__ == '__main__':
    main()
