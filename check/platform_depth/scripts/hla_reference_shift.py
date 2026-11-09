#!/usr/bin/env python3
"""Are the cases shifted toward ToMMo at HLA? A test of regional structure that the 10 PCs might miss (user
request 2026-10-30; HLA only: the one place where the control sources and the two Japanese panels differ as much
as cases and controls do).

If a share b of the cases came from where ToMMo 61KJPN was sampled and the rest resembled the controls, the
case-control difference of every HLA unit would follow the ToMMo-control difference with slope b; without such a
shift, b = 0. Units are ToMMo's reporting units (P groups or alleles; hla.typing.v2 05.exports/
allele_pgroup_map.tsv, in_reference = 1). The slope is fitted on class I and DP (A, B, C, E, F, G, DPA1, DPB1;
p >= 0.01 in controls or ToMMo); DR-DQ carries the DQB1 57 signal and is predicted from the fit. Three versions:

  frequency           d_case = p_case - p_ctrl on d_ref = p_ToMMo - p_ctrl (our two-field dosages summed into the
                      units; the published pooled random-model samples). Both differences contain p_ctrl, so each
                      of 100 splits takes d_case from one random half of the controls and d_ref from the other.
                      Weights 1/(p(1-p)) at the controls' and ToMMo's mean. P: permutation of case labels
                      (one-sided, b > 0); CI: bootstrap of cases and controls (10 splits per replicate).
  adjusted            the published random-model log OR of each tested allele (sex + 10 PCs; 06.retest/original,
                      equal to assoc_hla.v2 by gate G4) on delta = logit p_ToMMo - logit p_ctrl, weights 1/SE^2,
                      the shared control sampling variance 1/(2 n p (1-p)) removed from numerator and denominator.
  adjusted, given 57V the same with each allele's log OR given DQB1:57:V (assoc_hla.v2's conditional scan), so the
                      alleles' linkage to the signal does not enter the slope. With a real 57 V effect the
                      unconditional slope overstates the shift; with a regional 57 V difference the conditional one
                      understates it: the truth lies between them.
  CIs of the adjusted slopes: bootstrap over units (2,000); P: share of replicates <= 0.
For DQB1:57:V: predicted = b x d (frequency) or b x delta (log-odds); share explained = predicted / observed; the
residual (observed - predicted) / SE and its P.

Outputs (<out-prefix>.*):
  fit.tsv      per cohort x version: units, b, 95% CI, P, and DQB1:57:V observed, predicted, share, residual z, P
  groups.tsv   per cohort x unit: frequencies, d_case, d_ref, in the fit, predicted d_case (frequency version)
  signal.tsv   per cohort: DQB1:57:V and its P groups on the frequency scale (observed vs predicted)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import logit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402

FIT_GENES = ['A', 'B', 'C', 'E', 'F', 'G', 'DPA1', 'DPB1']
DRDQ = ['DRB1', 'DQA1', 'DQB1']


def residue57(exports):
    """DQB1 reporting unit -> residue at IMGT position 57 (one residue per P group; 'mixed' otherwise)."""
    ref = pd.read_csv(exports / 'residue_reference' / 'DQB1.tsv', sep='\t', dtype=str, usecols=['allele', '57'])
    mm = pd.read_csv(exports / 'allele_match_map.tsv', sep='\t', dtype=str)
    mm = mm[mm.gene == 'DQB1'].merge(ref.rename(columns={'allele': 'matched_to'}), on='matched_to', how='left')
    pg = pd.read_csv(exports / 'allele_pgroup_map.tsv', sep='\t', dtype=str)
    x = pg[pg.gene == 'DQB1'].merge(mm[['allele', '57']], on='allele', how='left')
    return x.groupby('p_group')['57'].agg(lambda v: v.dropna().iloc[0] if v.dropna().nunique() == 1 else 'mixed')


def freqs(G, called, rows):
    num = np.nansum(G[rows], axis=0)
    den = 2 * called[rows].sum(axis=0)
    return np.divide(num, den, out=np.full(G.shape[1], np.nan), where=den > 0)


def slope_freq(G, called, case_rows, ctrl_rows, p_ref, w, fit, rng, splits):
    p_case = freqs(G, called, case_rows)
    bs = []
    for _ in range(splits):
        perm = rng.permutation(ctrl_rows)
        a, b = perm[: len(perm) // 2], perm[len(perm) // 2:]
        dc = p_case - freqs(G, called, a)
        dr = p_ref - freqs(G, called, b)
        ok = fit & np.isfinite(dc) & np.isfinite(dr)
        bs.append(np.sum(w[ok] * dc[ok] * dr[ok]) / np.sum(w[ok] * dr[ok] ** 2))
    return float(np.mean(bs))


def slope_adj(f, col):
    w = 1 / f.se ** 2
    return float((np.sum(w * f[col] * f.delta) - np.sum(w * f.v)) / (np.sum(w * f.delta ** 2) - np.sum(w * f.v)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exports', default=str(K.HLA_V2 / 'results' / '05.exports'))
    ap.add_argument('--conditional', default=str(K.ASSOC_HLA / 'results' / 'all'),
                    help="assoc_hla.v2 results/all (<cohort>/04.conditional/random/hla_given_signal.tsv)")
    ap.add_argument('--results', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--splits', type=int, default=100)
    ap.add_argument('--perms', type=int, default=1000)
    ap.add_argument('--boot', type=int, default=500)
    ap.add_argument('--boot-units', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20261030)
    a = ap.parse_args()
    ex, res = Path(a.exports), Path(a.results)

    pg = pd.read_csv(ex / 'allele_pgroup_map.tsv', sep='\t', dtype=str)
    pg = pg[(pg.in_reference == '1') & pg.gene.isin(FIT_GENES + DRDQ)]
    units = pg.drop_duplicates('p_group')[['gene', 'p_group', 'freq_reference']].reset_index(drop=True)
    units['freq_reference'] = units.freq_reference.astype(float)
    units['residue57'] = units.p_group.map(residue57(ex)).where(units.gene == 'DQB1')
    col = {u: i for i, u in enumerate(units.p_group)}
    v_units = units.index[(units.gene == 'DQB1') & (units.residue57 == 'V')]
    p_ref_v57 = float(units.freq_reference[v_units].sum())
    dos = pd.read_csv(ex / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    alle = [x for x in pg.allele if x in dos.columns]
    amap = pg.set_index('allele').p_group

    fits, groups, signal = [], [], []
    for c in K.COHORTS:
        rng = np.random.default_rng(a.seed)
        # ── frequency version ───────────────────────────────────────────────
        ph = pd.read_csv(res / '08.subsets' / c / 'subsets' / 'pooled.random.tsv', sep='\t', dtype=str)
        ids = [s for s in ph.IID if s in dos.index]
        is_case = ph.set_index('IID').loc[ids].PHENO1.eq('1').values
        D = dos.loc[ids, alle]
        G = np.zeros((len(ids), len(units)))
        called = np.zeros((len(ids), len(units)), dtype=bool)
        for gene, gidx in units.groupby('gene').groups.items():
            ga = [x for x in alle if x.split('*')[0] == gene]
            ok = D[ga].notna().all(axis=1).values
            for x in ga:
                G[:, col[amap[x]]] += D[x].fillna(0).values
            for j in gidx:
                called[:, j] = ok
        G[~called] = np.nan
        case_rows, ctrl_rows = np.where(is_case)[0], np.where(~is_case)[0]
        p_ref = units.freq_reference.values
        p_ctrl, p_case = freqs(G, called, ctrl_rows), freqs(G, called, case_rows)
        pbar = (p_ctrl + p_ref) / 2
        w = 1 / (pbar * (1 - pbar))
        fit = units.gene.isin(FIT_GENES).values & (np.fmax(p_ctrl, p_ref) >= 0.01)
        b = slope_freq(G, called, case_rows, ctrl_rows, p_ref, w, fit, rng, a.splits)
        sub = max(1, a.splits // 10)
        null = [slope_freq(G, called, *(lambda l: (np.where(l)[0], np.where(~l)[0]))(rng.permutation(is_case)),
                           p_ref, w, fit, rng, sub) for _ in range(a.perms)]
        boot = [slope_freq(G, called, rng.choice(case_rows, len(case_rows)), rng.choice(ctrl_rows, len(ctrl_rows)),
                           p_ref, w, fit, rng, sub) for _ in range(a.boot)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        n_case, n_ctrl = len(case_rows), len(ctrl_rows)
        pc, pk = float(np.nansum(p_case[v_units])), float(np.nansum(p_ctrl[v_units]))
        dc, dr = pc - pk, p_ref_v57 - pk
        se = np.sqrt(pc * (1 - pc) / (2 * n_case) + pk * (1 - pk) / (2 * n_ctrl))
        fits.append({'cohort': c, 'version': 'frequency', 'units': int(fit.sum()), 'b': b, 'b_l95': lo, 'b_u95': hi,
                     'p_b_gt_0': (1 + np.sum(np.array(null) >= b)) / (1 + len(null)),
                     'v57_observed': dc, 'v57_predicted': b * dr, 'v57_predicted_at_b_u95': hi * dr,
                     'v57_share_explained': b * dr / dc, 'v57_residual_z': (dc - b * dr) / se,
                     'v57_residual_p': 2 * stats.norm.sf(abs((dc - b * dr) / se))})
        groups.append(units.assign(cohort=c, p_case=p_case, p_ctrl=p_ctrl, d_case=p_case - p_ctrl,
                                   d_ref=p_ref - p_ctrl, in_fit=fit, predicted_d_case=b * (p_ref - p_ctrl)))
        for name, idx in [(RL.V57, list(v_units))] + [(units.p_group[j], [j]) for j in v_units
                                                           if np.fmax(p_ctrl[j], p_ref[j]) >= 0.01]:
            q_c, q_k, q_r = np.nansum(p_case[idx]), np.nansum(p_ctrl[idx]), np.nansum(p_ref[idx])
            s_ = np.sqrt(q_c * (1 - q_c) / (2 * n_case) + q_k * (1 - q_k) / (2 * n_ctrl))
            signal.append({'cohort': c, 'unit': name, 'p_case': q_c, 'p_ctrl': q_k, 'p_tommo': q_r,
                           'd_case': q_c - q_k, 'd_ref': q_r - q_k, 'predicted_d_case': b * (q_r - q_k),
                           'share_explained': b * (q_r - q_k) / (q_c - q_k) if q_c != q_k else np.nan,
                           'z_residual': (q_c - q_k - b * (q_r - q_k)) / s_})

        # ── adjusted versions (published log ORs) ───────────────────────────
        t = pd.read_csv(res / '06.retest' / 'original' / c / 'assoc.random.allele.tsv', sep='\t')
        cond = pd.read_csv(Path(a.conditional) / c / '04.conditional' / 'random' / 'hla_given_signal.tsv', sep='\t')
        cond = cond[(cond['class'] == 'allele') & (cond.conditioned_on == RL.V57)].drop_duplicates('marker')
        cond = cond.set_index('marker')
        pgi = pg.set_index('allele')
        t = t[t.marker.isin(pgi.index)].copy()
        t['gene'] = t.marker.str.split('*').str[0]
        t['unit'] = pgi.loc[t.marker, 'p_group'].values
        t['p_ref'] = pgi.loc[t.marker, 'freq_reference'].astype(float).values
        t = t[~t.unit.duplicated(keep=False)]           # one tested allele per ToMMo unit
        t = t[np.isfinite(t.beta) & np.isfinite(t.se) & (t.p_ref > 0) & (t.freq_control > 0)].copy()
        t['beta_given'] = np.log(cond.reindex(t.marker).or_given.astype(float).values)
        t['delta'] = logit(t.p_ref) - logit(t.freq_control)
        t['v'] = 1 / (2 * t.called_control.astype(float) * t.freq_control * (1 - t.freq_control))
        f = t[t.gene.isin(FIT_GENES) & (np.fmax(t.freq_control, t.p_ref) >= 0.01)]
        rr = pd.read_csv(res / '06.retest' / 'original' / c / 'assoc.random.residue.tsv', sep='\t').set_index('marker')
        v57 = rr.loc[RL.V57]
        d57 = logit(p_ref_v57) - logit(v57.freq_control)
        for version, colname in [('adjusted', 'beta'), ('adjusted, given 57V', 'beta_given')]:
            ff = f.dropna(subset=[colname])
            bb = slope_adj(ff, colname)
            bs = np.array([slope_adj(ff.sample(len(ff), replace=True, random_state=rng.integers(1 << 31)), colname)
                           for _ in range(a.boot_units)])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            z = (v57.beta - bb * d57) / v57.se
            fits.append({'cohort': c, 'version': version, 'units': len(ff), 'b': bb, 'b_l95': lo, 'b_u95': hi,
                         'p_b_gt_0': (1 + np.sum(bs <= 0)) / (1 + len(bs)), 'v57_observed': v57.beta,
                         'v57_predicted': bb * d57, 'v57_predicted_at_b_u95': hi * d57,
                         'v57_share_explained': bb * d57 / v57.beta, 'v57_residual_z': z,
                         'v57_residual_p': 2 * stats.norm.sf(abs(z)),
                         'v57_residual_p_at_b_u95': 2 * stats.norm.sf(abs((v57.beta - hi * d57) / v57.se))})
        for r in [x for x in fits if x['cohort'] == c]:
            print(f"[hla_reference_shift] {c} {r['version']}: b {r['b']:.3f} ({r['b_l95']:.3f} to {r['b_u95']:.3f}), "
                  f"P {r['p_b_gt_0']:.3g}; 57V share {r['v57_share_explained']:.2f}, residual P {r['v57_residual_p']:.2g}")

    pd.DataFrame(fits).to_csv(f'{a.out_prefix}.fit.tsv', sep='\t', index=False)
    pd.concat(groups).to_csv(f'{a.out_prefix}.groups.tsv', sep='\t', index=False)
    pd.DataFrame(signal).to_csv(f'{a.out_prefix}.signal.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
