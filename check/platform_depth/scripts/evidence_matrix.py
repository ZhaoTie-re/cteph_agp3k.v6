#!/usr/bin/env python3
"""Evidence matrix: every signal x every check (x cohort), with its metric, its rule and its reading (report stage).

No signal is removed here. Each check reads 'consistent', 'flag', 'not testable' or 'not applicable' (the array
check: 'consistent', 'not checkable', or 'flag' when contradicted). Every check is read in each of the
three cohorts and flags a signal when it flags in any one of them, except the array check of the gene sets: one
'reported' cell for all three cohorts (user 2026-11-05: the cohorts are nested subsets), the case carriers of the
drivers confirmed on the imputed arrays, counted in full mainland, which contains the other two; a description, no rule.
The gene sets read are the M/H sets (checks.GENE_SETS_READ, the main analysis); the L/M/H re-tests stay in the stage
tables as data and are not read here. The checks,
their rules and thresholds are the registry in checks.py (proposed, not yet confirmed by the user); every number
comes from a published table of this component.

  evidence.long.tsv     signal x check x cohort: metric, value, reading
  evidence.loo.tsv      every leave-one-out estimate behind 'lopo' and 'sources_loo': z observed and expected, the
                        pooled z, the difference test P, and the calibration of the rule (below)
  evidence.loo_calibration.tsv  per check x cohort: testable readings, flags observed, flags expected by chance

Calibration of the leave-one-out rule (a diagnostic, not a rule): the estimate without a group is nested in the
pooled one, so if the effect is the same in every group, z without the group, given the pooled estimate, is
normal with mean z expected and SD sqrt(1 - (z expected / z pooled)^2) (the winner's curse inflates both and does
not bias the comparison). The chance of a flag of a testable reading follows from it; their sum is the number of
flags expected with no artifact at all.
  evidence.signals.tsv  one row per signal: discovery, the reading of every check over the three cohorts, the flags
  evidence.rules.tsv    the rule of every check
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import checks as CK  # noqa: E402
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402

zp = RL.z_from_p


class Matrix:
    def __init__(self):
        self.rows, self.loo = [], []

    def add_loo(self, sig, kind, check, cohort, rows, z_pooled=np.nan, p_diff=None):
        for i, (lab, zo, ze) in enumerate(rows):
            sd = np.sqrt(max(1 - (ze / z_pooled) ** 2, 0)) if np.isfinite(z_pooled) and z_pooled != 0 else np.nan
            testable = np.isfinite(ze) and abs(ze) >= CK.LOO_Z_EXPECTED
            chance = (norm.cdf((CK.LOO_Z_OBSERVED - abs(ze)) / sd) if testable and sd > 0 else np.nan)
            self.loo.append({'signal': sig, 'kind': kind, 'check': check, 'cohort': cohort, 'left_out': lab,
                             'z_observed': zo, 'z_expected': ze, 'z_pooled': z_pooled, 'sd_conditional': sd,
                             'u': (zo - ze) / sd if sd and sd > 0 else np.nan,
                             'p_diff': p_diff[i] if p_diff is not None else np.nan,
                             'testable': bool(testable), 'chance_of_flag': chance})
        return RL.loo_reading(rows)

    def add(self, sig, kind, check, cohort, reading, metric, value=np.nan):
        self.rows.append({'signal': sig, 'kind': kind, 'check': check, 'cohort': cohort, 'reading': reading,
                          'metric': metric, 'value': value})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', type=Path, default=K.ROOT / 'results')
    ap.add_argument('--label', default='downsampling')
    ap.add_argument('--ds-summary', type=Path, help='folder of ds.* (default: results/07.downsampling/<label>/summary)')
    ap.add_argument('--array-tables', type=Path, help='folder of array.* (default: results/11.array)')
    ap.add_argument('--control-tables', type=Path, help='folder of controls.* (default: results/12.controls)')
    ap.add_argument('--followup-tables', type=Path, help='unused: the follow-up counts are descriptive, not a check')
    ap.add_argument('--out-dir', type=Path, default=Path('.'))
    a = ap.parse_args()
    R, lab = a.results, a.label
    S = a.ds_summary or R / '07.downsampling' / lab / 'summary'
    AR = a.array_tables or R / '11.array'
    CT = a.control_tables or R / '12.controls'
    M = Matrix()

    sig = RL.read_signals(R)
    saige = RL.saige_signals(sig)
    n_saige = len(saige)
    gsets = [gs for gs in RL.gene_sets(sig) if gs[1] in CK.GENE_SETS_READ]   # M/H only; L/M/H stays data
    q_threshold = CK.Q_ALPHA / n_saige

    dsg = pd.read_csv(S / 'ds.gwas.tsv', sep='\t')
    dsg = dsg[(dsg.role == 'Lead') & (dsg.imputation == 'best_guess')]
    dsh = pd.read_csv(S / 'ds.hla_dqb1_57v.tsv', sep='\t')
    dsp = pd.read_csv(S / 'ds.positions.tsv', sep='\t')
    dsn = pd.read_csv(S / 'ds.genes.tsv', sep='\t')
    REPS = RL.replicate_arms(dsg.columns, 'p_')
    for t, prefix in ((dsh, 'p_'), (dsp, 'p_'), (dsn, 'pvalue_')):
        if RL.replicate_arms(t.columns, prefix) != REPS:
            sys.exit('ABORT: the ds.* tables hold different replicate arms')
    n_reps = len(REPS)
    arr = pd.read_csv(AR / 'array.signals.tsv', sep='\t')
    chet = pd.read_csv(CT / 'controls.heterogeneity.tsv', sep='\t', dtype={'variant_id': str})
    cloo = pd.read_csv(CT / 'controls.loo.tsv', sep='\t', dtype={'variant_id': str})
    cgen = pd.read_csv(CT / 'controls.genes.tsv', sep='\t')
    cpos = pd.read_csv(CT / 'controls.positions.tsv', sep='\t', dtype=str)

    for c in K.COHORTS:
        het = pd.read_csv(RL.platform_table(R, c, 'heterogeneity'), sep='\t', dtype={'variant_id': str}).set_index('variant_id')
        lo = pd.read_csv(RL.platform_table(R, c, 'lopo'), sep='\t', dtype={'variant_id': str})
        ch = chet[chet.cohort == c].set_index('variant_id')
        cl = cloo[cloo.cohort == c]
        ar = arr[arr.cohort == c].set_index('signal')

        # ── SAIGE signals: the GWAS leads and DQB1:57:V ──────────────────────────────
        for unit, vid in saige:
            kind = 'hla_residue' if unit == RL.V57 else 'gwas'
            if kind == 'gwas':
                d = dsg[(dsg.unit == unit) & (dsg.cohort == c)].iloc[0]
                po = d.p_original
                lvl = 'genome-wide' if po < 5e-8 else 'suggestive' if po < 1e-5 else 'below 1e-5'
                M.add(unit, kind, 'discovery', c, lvl, f'P {po:.2e}', po)
                pr = np.median([d[f'p_{r}'] for r in REPS])
                nsig = int(d['reps_p_lt_1e-5'])
                extra = f'; replicates P < 1e-5: {nsig}/{n_reps}, < 5e-8: {int(d["reps_p_lt_5e-8"])}/{n_reps}'
            else:
                d = dsh[dsh.cohort == c].iloc[0]
                po = d.p_original
                gec = po / d.p_over_gec_original
                M.add(unit, kind, 'discovery', c, 'significant' if po < gec else 'not significant',
                      f'P {po:.2e} (GEC {gec:.1e})', po)
                pr = np.median([d[f'p_{r}'] for r in REPS])
                nsig = int(sum(d[f'p_over_gec_{r}'] < 1 for r in REPS))
                extra = f'; replicates significant at GEC: {nsig}/{n_reps}'
            dz = float(zp(po) - zp(pr))
            M.add(unit, kind, 'depth', c, 'flag' if dz > CK.DEPTH_DZ else 'consistent',
                  f'|z| original {zp(po):.2f}, replicate median {zp(pr):.2f}{extra}', dz)

            h = het.loc[vid]
            maf = min(h.af_pooled, 1 - h.af_pooled)
            if maf < CK.Q_MIN_MAF:
                M.add(unit, kind, 'platform_Q', c, 'not applicable', f'MAF {maf:.3f}; Q P {h.p_Q_empirical:.3f}, I2 {h.I2:.2f}',
                      h.p_Q_empirical)
            else:
                M.add(unit, kind, 'platform_Q', c, 'flag' if h.p_Q_empirical < q_threshold else 'consistent',
                      f'Q P {h.p_Q_empirical:.3f} (threshold {q_threshold:.4f}), I2 {h.I2:.2f}', h.p_Q_empirical)
            re_txt = (f'OR pooled {np.exp(h.beta_pooled):.2f}, random effects {np.exp(h.beta_re):.2f} '
                      f'({np.exp(h.re_l95):.2f}-{np.exp(h.re_u95):.2f})')
            if maf < CK.Q_MIN_MAF:          # as Q: SAIGE's one-step stratum beta is unreliable in small strata
                M.add(unit, kind, 'platform_RE', c, 'not applicable', f'MAF {maf:.3f}; {re_txt}', np.exp(h.beta_re))
            else:
                M.add(unit, kind, 'platform_RE', c, 'flag' if np.sign(h.beta_re) != np.sign(h.beta_pooled) else 'consistent',
                      re_txt, np.exp(h.beta_re))
            x_lo = lo[lo.variant_id == vid]
            r, txt = M.add_loo(unit, kind, 'lopo', c, [(x.left_out, x.z_obs, x.z_expected) for x in x_lo.itertuples()],
                               h.beta_pooled / h.se_pooled, list(x_lo.p_diff))
            M.add(unit, kind, 'lopo', c, r, txt)

            if vid in ch.index:
                q = ch.loc[vid]
                if maf < CK.Q_MIN_MAF:
                    M.add(unit, kind, 'sources_Q', c, 'not applicable', f'MAF {maf:.3f}; Q P {q.p_Q:.3f}, I2 {q.I2:.2f}', q.p_Q)
                else:
                    M.add(unit, kind, 'sources_Q', c, 'flag' if q.p_Q < q_threshold else 'consistent',
                          f'Q P {q.p_Q:.3f}, I2 {q.I2:.2f}', q.p_Q)
            x_cl = cl[cl.variant_id == vid]
            r, txt = M.add_loo(unit, kind, 'sources_loo', c, [(x.left_out, x.z_obs, x.z_expected) for x in x_cl.itertuples()],
                               ch.loc[vid].beta_pooled / ch.loc[vid].se_pooled if vid in ch.index else np.nan,
                               list(x_cl.p_diff))
            M.add(unit, kind, 'sources_loo', c, r, txt)

            if unit in ar.index:
                v = ar.loc[unit]
                rd = 'flag' if 'contradict' in str(v.array_verdict) else CK.ARRAY_READING.get(v.array_verdict, v.array_verdict)
                M.add(unit, kind, 'array', c, rd, f'{v.array_verdict} ({v.strongest_evidence})')

        # ── the residue position DQB1 57 (fixed omnibus) ─────────────────────────────
        pos = 'DQB1 position 57'
        d = dsp[(dsp.unit == RL.POSITION57) & (dsp.cohort == c)].iloc[0]
        po, thr = d.p_original, d.p_original / d.p_over_threshold_original
        M.add(pos, 'hla_position', 'discovery', c, 'significant' if po < thr else 'not significant',
              f'P {po:.2e} (Bonferroni {thr:.1e})', po)
        pr = np.median([d[f'p_{r}'] for r in REPS])
        dz = float(zp(po) - zp(pr))
        M.add(pos, 'hla_position', 'depth', c, 'flag' if dz > CK.DEPTH_DZ else 'consistent',
              f'|z| original {zp(po):.2f}, replicate median {zp(pr):.2f}; replicates significant: '
              f'{int(d.reps_bonferroni)}/{n_reps}', dz)
        cp = cpos[cpos.cohort == c] if 'cohort' in cpos.columns else cpos.iloc[0:0]
        cp = cp[cp.unit == RL.POSITION57] if len(cp) else cp
        if len(cp):
            ps = {k[2:]: float(v) for k, v in cp.iloc[0].items() if k.startswith('p_ctrlwithout_')}
            worst = max(ps, key=ps.get)
            M.add(pos, 'hla_position', 'sources_loo', c, 'reported',
                  f'largest P without one source: -{worst.replace("ctrlwithout_", "")} {ps[worst]:.1e}', ps[worst])
        ps = {}
        for p in RL.CASE_SLUGS:
            f = RL.omnibus_file(R, f'without_{p}')
            if f.exists():
                o = pd.read_csv(f, sep='\t', dtype=str)
                o = o[(o.cohort == c) & (o.unit == RL.POSITION57) & (o.case_set == 'all')]
                if len(o):
                    ps[p] = float(o.p.iloc[0])
        if ps:
            worst = max(ps, key=ps.get)
            M.add(pos, 'hla_position', 'lopo', c, 'reported',
                  f'largest P without one platform: -{worst} {ps[worst]:.1e} (Bonferroni {thr:.1e})', ps[worst])

        # ── genes: Tier 1-2, each tiered set ──────────────────────────────────────────
        for g, st in gsets:
            name = f'{g} ({st})'
            x = dsn[(dsn.gene == g) & (dsn.stratum == st) & (dsn.cohort == c)]
            bon = x.bonferroni_threshold.iloc[0]
            best = x.loc[x.pvalue_original.idxmin()]
            M.add(name, 'gene', 'discovery', c, 'significant' if best.pvalue_original < bon else 'not significant',
                  '; '.join(f'{r.method} P {r.pvalue_original:.1e}' for r in x.itertuples()) + f' (Bonferroni {bon:.1e})',
                  best.pvalue_original)
            dzs = []
            for r in x.itertuples():
                pr = np.median([getattr(r, f'pvalue_{k}') for k in REPS])
                dzs.append((r.method, float(zp(r.pvalue_original) - zp(pr)), int(r.reps_bonferroni),
                            int(r.carriers_case_rep_min), int(r.carriers_case_rep_max), int(r.carriers_case_original)))
            dzmax = max(v[1] for v in dzs)
            M.add(name, 'gene', 'depth', c, 'flag' if dzmax > CK.DEPTH_DZ else 'consistent',
                  '; '.join(f'{m}: |z| change {v:+.2f}, Bonferroni {b}/{n_reps}' for m, v, b, *_ in dzs) +
                  f'; case carriers {dzs[0][5]} -> {dzs[0][3]}-{dzs[0][4]}', dzmax)
            gp = pd.read_csv(RL.gene_subset_file(R, c, 'pooled', st), sep='\t')
            p0 = gp[(gp.gene == g) & (gp.method == 'cmc')].iloc[0]
            rows = []
            for p in RL.CASE_SLUGS:
                w = pd.read_csv(RL.gene_subset_file(R, c, f'without_{p}', st), sep='\t')
                w = w[(w.gene == g) & (w.method == 'cmc')]
                if len(w) and w.se.iloc[0] > 0:
                    rows.append((p, w.beta.iloc[0] / w.se.iloc[0], p0.beta / w.se.iloc[0]))
            r, txt = M.add_loo(name, 'gene', 'lopo', c, rows, p0.beta / p0.se if p0.se > 0 else np.nan)
            M.add(name, 'gene', 'lopo', c, r, txt)
            for chk in ('platform_Q', 'platform_RE', 'sources_Q'):
                M.add(name, 'gene', chk, c, 'not applicable', 'carriers per stratum too few for a per-stratum estimate')
            cg = cgen[(cgen.cohort == c) & (cgen.gene == g) & (cgen.stratum == st) & (cgen.method == 'cmc')].set_index('subset')
            rows = []
            for s, k in zip(RL.SOURCES, RL.SOURCE_WITHOUT):
                if k in cg.index and cg.loc[k, 'se'] > 0:
                    rows.append((s, cg.loc[k, 'beta'] / cg.loc[k, 'se'], cg.loc['pooled', 'beta'] / cg.loc[k, 'se']))
            zp0 = cg.loc['pooled', 'beta'] / cg.loc['pooled', 'se'] if 'pooled' in cg.index and cg.loc['pooled', 'se'] > 0 else np.nan
            r, txt = M.add_loo(name, 'gene', 'sources_loo', c, rows, zp0)
            M.add(name, 'gene', 'sources_loo', c, r, txt)
            if c == CK.ARRAY_GENE_COHORT and g in ar.index:    # one cell for the three cohorts (nested subsets)
                v = ar.loc[g]
                rd = 'flag' if 'contradict' in str(v.array_verdict) else CK.ARRAY_READING.get(v.array_verdict, v.array_verdict)
                M.add(name, 'gene', 'array', c, rd,
                      (f'case carriers confirmed on the imputed arrays: {v.case_carriers_recovered.replace("cases ", "")}'
                       if v.array_verdict == 'carriers counted' else 'no driver in the imputed arrays') +
                      ' (full mainland, which contains the narrow and intermediate cohorts)')

    long = pd.DataFrame(M.rows)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(a.out_dir / 'evidence.long.tsv', sep='\t', index=False)
    loo = pd.DataFrame(M.loo)
    loo.to_csv(a.out_dir / 'evidence.loo.tsv', sep='\t', index=False)
    t = loo[loo.testable]
    flagged = RL.loo_flags(t.z_observed, t.z_expected)
    cal = (t.assign(flag=flagged).groupby(['check', 'cohort'])
           .agg(testable=('flag', 'size'), flags_observed=('flag', 'sum'), flags_expected=('chance_of_flag', 'sum'))
           .reset_index())
    cal.to_csv(a.out_dir / 'evidence.loo_calibration.tsv', sep='\t', index=False)
    pd.DataFrame([{'check': k, 'rule': v} for k, v in CK.RULES.items()]).to_csv(a.out_dir / 'evidence.rules.tsv', sep='\t', index=False)

    order = list(CK.RULES)
    out = []
    for (s, kind), x in long.groupby(['signal', 'kind'], sort=False):
        row = {'signal': s, 'kind': kind}
        flags = []
        for chk in order:
            y = x[x.check == chk]
            if not len(y):
                row[chk] = ''
                continue
            by = dict(zip(y.cohort, y.reading))
            if chk == 'discovery':
                row[chk] = ' / '.join(by.get(c, '-') for c in K.COHORTS)
                continue
            rd = CK.summary_reading([by[c] for c in K.COHORTS if c in by], chk)
            row[chk] = rd
            if rd == 'flag':
                flags.append(chk)
        row['flags'] = ', '.join(flags) if flags else 'none'
        out.append(row)
    sigs = pd.DataFrame(out)
    sigs.to_csv(a.out_dir / 'evidence.signals.tsv', sep='\t', index=False)
    print(f'[evidence_matrix] {len(sigs)} signals x {len(order)} checks; flagged: ' +
          ', '.join(f'{r.signal} ({r.flags})' for r in sigs.itertuples() if r.flags != 'none'))


if __name__ == '__main__':
    main()
