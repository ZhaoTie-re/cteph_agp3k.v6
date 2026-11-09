#!/usr/bin/env python3
"""Down-sampling: every signal in every cohort, ORIGINAL (published genotypes / typing) vs the f = 1
baseline vs the R replicates, with each test's pre-specified threshold.

  gwas       the locus leads and alternative leads (and every credible-set member in the long table),
             SAIGE best_guess (as published) and mean imputation: beta, OR, P; replicates as median and
             range; replicates reaching P < 1e-5 (suggestive) and < 5e-8 (genome-wide)
  genes      the Tier 1-2 genes x variant set x CMC / SKAT-O: P, OR, carriers by status; Bonferroni and
             BH (the gene's P substituted into its published scan)
  hla        every tested allele and residue (random model): OR, P, P / GEC threshold; n/3 cohorts
             significant at GEC per arm; DQB1:57:V itemised
  positions  the residue-position omnibus: chi2, P, P / Bonferroni threshold
  exposure   per cohort, how many of each signal's case carriers / case copies are among the down-sampled
             cases (all of them, the deep NovaSeq and DNBSeq-G400RS 30x cases, the DNBSeq-T7 cases): only
             those can move
  qc_reasons per cohort, the share of genotypes failing each QC rule (median over each group's genomes; the
             re-called arms also per platform of the down-sampled cases)
The down-sampled cases are deep_cases.tsv of the arm plan (role 'deep' there): every case platform but
DNBseq-G400RS 15x (user, 2026-11-03).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402

ASSOC_HLA = K.ASSOC_HLA / 'results' / 'all'


def spread(df, key, col, reps):
    """median / min / max of `col` over the replicate arms, per key."""
    r = df[df.arm.isin(reps)].groupby(key)[col]
    return pd.DataFrame({f'{col}_rep_median': r.median(), f'{col}_rep_min': r.min(), f'{col}_rep_max': r.max()})


RETEST = {}       # 'original': <retest>/original, 'label': <retest>/<run label> (set in main)


def load(pattern):
    """<retest>/original/<cohort>/<file> and <retest>/<run label>/<arm>/<cohort>/<file>."""
    rows = []
    files = [(('original',) + f.parts[-2:-1], f) for f in sorted(Path(RETEST['original']).glob(pattern))]
    files += [((f.parts[-3], f.parts[-2]), f) for f in sorted(Path(RETEST['label']).glob('*/' + pattern))]
    for (arm, cohort), f in files:
        d = pd.read_csv(f, sep='\t', dtype={'variant_id': str, 'marker': str})
        d['arm'], d['cohort'] = arm, cohort
        rows.append(d)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--retest', required=True, help='results/06.retest')
    ap.add_argument('--label', required=True, help='the run label of the arms (pilot / downsampling)')
    ap.add_argument('--signals', required=True)
    ap.add_argument('--deep', required=True)
    ap.add_argument('--samples', required=True)
    ap.add_argument('--orig-geno', required=True, help='results/04.genotypes/original (one folder per cohort)')
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    RETEST['original'] = Path(a.retest) / 'original'
    RETEST['label'] = Path(a.retest) / a.label
    dsc = pd.read_csv(a.deep, sep='\t', dtype=str)
    deep = set(dsc.sample_id)                                    # every down-sampled case
    sub = {'deep': set(dsc[dsc.platform.isin(K.DEEP_PLATFORMS)].sample_id),
           't7': set(dsc[dsc.platform == 'DNBSeq-T7 30x'].sample_id)}
    sig = pd.read_csv(a.signals, sep='\t', dtype=str)
    out = {}

    # ── GWAS ─────────────────────────────────────────────────────────────────
    g = load('*/gwas.*.tsv')
    reps = sorted(x for x in g.arm.unique() if x.startswith('r'))
    gs = sig[sig.kind == 'gwas'][['unit', 'unit_id', 'role', 'variant_id']].drop_duplicates('variant_id')
    g = g.merge(gs, on='variant_id', how='inner')
    g['mlog10p'] = -np.log10(g.p)
    key = ['unit', 'unit_id', 'role', 'variant_id', 'cohort', 'imputation']
    wide = g.pivot_table(index=key, columns='arm', values=['beta', 'p'], aggfunc='first')
    wide.columns = [f'{v}_{a}' for v, a in wide.columns]
    wide = wide.join(spread(g, key, 'beta', reps)).join(spread(g, key, 'mlog10p', reps))
    r = g[g.arm.isin(reps)].groupby(key)
    wide['reps'] = r.size()
    wide['reps_p_lt_1e-5'] = r.p.apply(lambda p: int((p < 1e-5).sum()))
    wide['reps_p_lt_5e-8'] = r.p.apply(lambda p: int((p < 5e-8).sum()))
    wide = wide.reset_index()
    wide.to_csv(f'{a.out_prefix}.gwas_long.tsv', sep='\t', index=False)
    wide[wide.role.isin(['Lead', 'Alt lead'])].to_csv(f'{a.out_prefix}.gwas.tsv', sep='\t', index=False)
    out['gwas_variants'] = int(wide.variant_id.nunique())

    # ── genes ────────────────────────────────────────────────────────────────
    gn = load('*/gene.*.tsv')
    if len(gn):
        gn['mlog10p'] = -np.log10(gn.pvalue)
        key = ['gene', 'stratum', 'method', 'cohort']
        w = gn.pivot_table(index=key, columns='arm', values=['pvalue', 'or', 'carriers_case', 'carriers_control'],
                           aggfunc='first')
        w.columns = [f'{v}_{a}' for v, a in w.columns]
        w = w.join(spread(gn, key, 'mlog10p', reps)).join(spread(gn, key, 'carriers_case', reps))
        rr = gn[gn.arm.isin(reps)].groupby(key)
        w['reps'] = rr.size()
        w['reps_bonferroni'] = rr.significant_bonferroni.sum()
        w['reps_bh'] = rr.significant_bh.sum()
        w['bonferroni_threshold'] = gn.groupby(key).bonferroni_threshold.first()
        w.reset_index().to_csv(f'{a.out_prefix}.genes.tsv', sep='\t', index=False)

    # ── HLA single markers and positions ─────────────────────────────────────
    h = load('*/assoc.random.*.tsv')
    if len(h):
        thr = {}
        for c in K.COHORTS:
            for k in ('allele', 'residue'):
                t = pd.read_csv(ASSOC_HLA / c / '01.markers' / k / f'effective_tests.{k}.tsv', sep='\t')
                thr[(c, k)] = float(t.p_gec.iloc[0])
        h['gec_threshold'] = [thr[(c, k)] for c, k in zip(h.cohort, h['class'])]
        h['p_over_gec'] = h.p / h.gec_threshold
        h['sig'] = h.sig_gec.astype(str).str.lower().eq('true')
        key = ['class', 'gene', 'marker', 'cohort']
        w = h.pivot_table(index=key, columns='arm', values=['or', 'p', 'p_over_gec'], aggfunc='first')
        w.columns = [f'{v}_{a}' for v, a in w.columns]
        hh = h.assign(mlog10p=-np.log10(h.p))
        w = w.join(spread(hh, key, 'mlog10p', reps)).join(spread(hh, key, 'p_over_gec', reps))
        rr = h[h.arm.isin(reps)].groupby(key)
        w['reps'] = rr.size()
        w['reps_gec'] = rr.sig.sum()
        w = w.reset_index()
        w.to_csv(f'{a.out_prefix}.hla_long.tsv', sep='\t', index=False)
        # n/3 cohorts significant at GEC, per arm, for every marker significant in any arm
        n3 = h.groupby(['arm', 'class', 'marker']).sig.sum().rename('n_cohorts_gec').reset_index()
        anysig = set(n3[n3.n_cohorts_gec > 0].marker)
        n3[n3.marker.isin(anysig)].pivot_table(index=['class', 'marker'], columns='arm', values='n_cohorts_gec',
                                               aggfunc='first').reset_index().to_csv(f'{a.out_prefix}.hla_n3.tsv',
                                                                                    sep='\t', index=False)
        w[w.marker == RL.V57].to_csv(f'{a.out_prefix}.hla_dqb1_57v.tsv', sep='\t', index=False)
    om = []
    ofiles = [('original', Path(RETEST['original']) / '_omnibus' / 'omnibus.residue.tsv')]
    ofiles += [(f.parts[-3], f) for f in sorted(Path(RETEST['label']).glob('*/_omnibus/omnibus.residue.tsv'))]
    for arm, f in ofiles:
        if not f.exists():
            continue
        d = pd.read_csv(f, sep='\t', dtype=str)
        d['arm'] = arm
        om.append(d)
    if om:
        om = pd.concat(om)
        om['p'] = om.p.astype(float)
        om['p_over_threshold'] = om.p / om.threshold.astype(float)
        om['mlog10p'] = -np.log10(om.p)
        key = ['unit', 'cohort']
        w = om.pivot_table(index=key, columns='arm', values=['p', 'p_over_threshold'], aggfunc='first')
        w.columns = [f'{v}_{a}' for v, a in w.columns]
        w = w.join(spread(om, key, 'mlog10p', reps))
        rr = om[om.arm.isin(reps)].groupby(key)
        w['reps'] = rr.size()
        w['reps_bonferroni'] = rr.sig_bonf.apply(lambda x: int(x.astype(str).str.lower().eq('true').sum()))
        w.reset_index().to_csv(f'{a.out_prefix}.positions.tsv', sep='\t', index=False)

    # ── exposure: case carriers / copies among the down-sampled cases (each cohort, original genotypes) ─────
    st = K.read_tsv(a.samples).set_index('sample_id')
    rows = []
    out['dqb1_57v_copies'] = {}
    for c in K.COHORTS:
        for name in ('gwas', 'gene'):
            fam, vids, dose = RL.alt_dosage(f'{a.orig_geno}/{c}/{name}')
            want = set(sig[sig.kind == name].variant_id)
            is_case = np.array([st.loc[s, 'status'] == 'case' if s in st.index else False for s in fam])
            masks = {'downsampled': np.array([s in deep for s in fam]),
                     **{k: np.array([s in v for s in fam]) for k, v in sub.items()}}
            for vid, d in zip(vids, dose):
                if vid not in want:
                    continue
                called = np.isfinite(d)
                alt_n = np.where(called, d, 0)
                r = {'cohort': c, 'kind': name, 'variant_id': vid,
                     'case_alt_copies': int(alt_n[is_case & called].sum()),
                     'case_carriers': int(((alt_n > 0) & is_case & called).sum())}
                for k, m in masks.items():
                    r[f'{k}_alt_copies'] = int(alt_n[m & called].sum())
                    r[f'{k}_carriers'] = int(((alt_n > 0) & m & called).sum())
                rows.append(r)
        hm = pd.read_csv(K.HLA_V2 / 'results' / '06.assoc_inputs' / c / 'residue_dosage.tsv', sep='\t',
                         dtype={'sample_id': str}).set_index('sample_id')
        v = hm[RL.V57]
        cases = [s for s in hm.index if st.loc[s, 'status'] == 'case']
        out['dqb1_57v_copies'][c] = {'cases': float(v.loc[cases].sum()),
                                     'downsampled_cases': float(v.loc[[s for s in cases if s in deep]].sum()),
                                     **{f'{k}_cases': float(v.loc[[s for s in cases if s in m]].sum())
                                        for k, m in sub.items()}}
    ex = pd.DataFrame(rows)
    ex['downsampled_share_of_case_alt_copies'] = ex.downsampled_alt_copies / ex.case_alt_copies.where(ex.case_alt_copies > 0)
    ex.to_csv(f'{a.out_prefix}.exposure.tsv', sep='\t', index=False)

    # ── genotypes failing a QC rule, per cohort: the median share over each group's genomes (as the
    #    QC_REASONS summary, which counts the full-mainland members) ────────────────────────────────────────
    qr = Path(a.retest).parent / '07.downsampling' / a.label / 'qc_reasons' / 'qc_reasons.samples.tsv'
    if qr.exists():
        q = pd.read_csv(qr, sep='\t', dtype={'sample_id': str})
        rules = ['LowGQ', 'LowDP', 'ABB_outlier', 'ABB_NaN', 'any']
        parts = []
        q['platform'] = q.sample_id.map(st['platform'])
        for c in K.COHORTS:
            x = q[q.sample_id.map(st[f'in_{c}']) == '1']
            m = x.groupby(['source', 'group', 'arm'])[rules].median()
            m['samples'] = x.groupby(['source', 'group', 'arm']).size()
            parts.append(m.reset_index().assign(cohort=c, platform='all'))
            y = x[x.source == 'recalled']                    # the re-called arms, per platform as well
            m = y.groupby(['source', 'group', 'platform', 'arm'])[rules].median()
            m['samples'] = y.groupby(['source', 'group', 'platform', 'arm']).size()
            parts.append(m.reset_index().assign(cohort=c))
        qq = pd.concat(parts, ignore_index=True)[['cohort', 'source', 'group', 'platform', 'arm'] + rules + ['samples']]
        qq.to_csv(f'{a.out_prefix}.qc_reasons.tsv', sep='\t', index=False)
    out['replicates'] = reps
    Path(f'{a.out_prefix}.summary.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
