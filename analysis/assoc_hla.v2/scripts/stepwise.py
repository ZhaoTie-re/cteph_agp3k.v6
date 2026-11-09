#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Independent signals within one marker class for one (case set,
#           model), selected JOINTLY over the three cohorts, with the scan's own
#           engine (engine.py) on each cohort's pooled PLINK set.
#
#           Step 0 is the marginal scan. At step k every class marker is tested
#           conditional on the leads entered so far (in each cohort, the leads
#           that cohort tested). A marker is a CANDIDATE if it is significant
#           at the cohort's LD-adjusted (GEC) threshold in at least one cohort.
#           The step's lead is the candidate whose LARGEST P over the three
#           cohorts is smallest (a marker a cohort did not test counts P = 1
#           there): the marker most consistently associated, not the one with
#           the single smallest P. It ENTERS as an independent signal; the walk
#           stops when no marker is a candidate, after at most --max-signals
#           signals.
#
#           A signal's MEMBERS are the markers significant (GEC, any cohort)
#           before its lead entered and in no cohort after: the markers it
#           explains, its lead included. Each signal is reported with its
#           cross-cohort consistency at entry (conditional on the earlier
#           signals) and marginally, under the three thresholds; for a residue
#           lead, the tested alleles of its gene that carry the residue
#           (hla.typing.v2 allele_match_map + residue_reference).
#
#           Outputs: steps.<class>.tsv (one row per step), stepwise.<class>.tsv
#           (every marker at every step in every cohort), signals.<class>.tsv
#           (one row per entered signal), leads.<class>.txt (the entered leads).
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402
from engine import Engine              # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--model', required=True, choices=V.MODELS)
    p.add_argument('--class', dest='cls', required=True, choices=V.CLASSES)
    p.add_argument('--cohorts', required=True, help='comma list; every list below in this order')
    p.add_argument('--bed', nargs='+', required=True, help='pooled.bed per cohort (bim/fam beside it)')
    p.add_argument('--markers', nargs='+', required=True, help='pooled_markers.tsv per cohort')
    p.add_argument('--coding', nargs='+', required=True, help='pheno_coding.json per cohort')
    p.add_argument('--pheno-fixed', nargs='+', required=True)
    p.add_argument('--keep', nargs='+', required=True)
    p.add_argument('--rda', nargs='+', required=True)
    p.add_argument('--vr', nargs='+', required=True)
    p.add_argument('--thresholds', nargs='+', required=True, help='effective_tests.<class>.tsv per cohort')
    p.add_argument('--allele-map', required=True, help='hla.typing.v2 allele_match_map.tsv')
    p.add_argument('--residue-ref', required=True, help='hla.typing.v2 residue_reference/ directory')
    p.add_argument('--covar-name', default='SEX,PC1_AVG-PC10_AVG')
    p.add_argument('--plink2', required=True)
    p.add_argument('--saige-rscript', required=True)
    p.add_argument('--saige-step2', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--max-signals', type=int, default=5)
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def consistency(p, orr, thr):
    """Cohorts with p < thr and the majority direction (see compare.py)."""
    sig = [pd.notna(x) and x < t for x, t in zip(p, thr)]
    up = sum(1 for s, o in zip(sig, orr) if s and o > 1)
    down = sum(1 for s, o in zip(sig, orr) if s and o < 1)
    return max(up, down)


class Carriers:
    """Tested alleles of a gene that carry a residue at an IMGT position."""

    def __init__(self, match_map, ref_dir, alleles):
        self.map = pd.read_csv(match_map, sep='\t', dtype=str).set_index('allele')['matched_to'].to_dict()
        self.ref_dir, self.alleles, self.ref = Path(ref_dir), alleles, {}

    def __call__(self, marker):
        gene, pos, aa = V.parse_marker(marker, 'residue')
        if gene not in self.ref:
            self.ref[gene] = pd.read_csv(self.ref_dir / f'{gene}.tsv', sep='\t', dtype=str).set_index('allele')
        ref = self.ref[gene]
        out = []
        for al in sorted((x for x in self.alleles if x.startswith(f'{gene}*')), key=lambda n: V.sort_key(n, 'allele')):
            full = self.map.get(al)
            if full in ref.index and str(pos) in ref.columns and ref.at[full, str(pos)] == aa:
                out.append(al)
        return out


def main():
    a = parse_args()
    out = Path(a.out_dir)
    cohorts = a.cohorts.split(',')
    lists = [a.bed, a.markers, a.coding, a.pheno_fixed, a.keep, a.rda, a.vr, a.thresholds]
    if any(len(x) != len(cohorts) for x in lists):
        sys.exit('ERROR: every per-cohort list must have one entry per cohort')
    short = [V.COHORT_SHORT[c] for c in cohorts]

    eng, members_of, thr, alleles = {}, {}, {}, set()
    for i, c in enumerate(cohorts):
        mk = pd.read_csv(a.markers[i], sep='\t', dtype={'marker': str})
        if set(mk['cohort']) != {c} or set(mk['case_set']) != {a.case_set}:
            sys.exit(f'ERROR: {a.markers[i]} is not {a.case_set} / {c}')
        th = pd.read_csv(a.thresholds[i], sep='\t').iloc[0]
        if (th['case_set'], th['cohort'], th['class']) != (a.case_set, c, a.cls):
            sys.exit(f'ERROR: {a.thresholds[i]} is not {a.case_set} / {c} / {a.cls}')
        thr[c] = {t: float(th[f'p_{t}']) for t in V.THRESHOLDS}
        members_of[c] = mk.loc[mk['class'] == a.cls, 'marker'].tolist()
        alleles |= set(mk.loc[mk['class'] == 'allele', 'marker'])
        cfg = SimpleNamespace(model=a.model, bfile=a.bed[i][:-4], keep=a.keep[i], pheno_fixed=a.pheno_fixed[i],
                              covar_name=a.covar_name, rda=a.rda[i], vr=a.vr[i], plink2=a.plink2,
                              saige_rscript=a.saige_rscript, saige_step2=a.saige_step2, threads=a.threads)
        want = json.load(open(a.coding[i]))[a.model]
        eng[c] = Engine(cfg, dict(zip(mk['marker'], mk['file_order'])), want, prefix=f'{short[i]}_')
    universe = sorted(set().union(*members_of.values()), key=lambda n: V.sort_key(n, a.cls))
    carriers = Carriers(a.allele_map, a.residue_ref, alleles) if a.cls == 'residue' else None
    key = {'case_set': a.case_set, 'model': a.model, 'class': a.cls}

    cond, prev_sig, marg = [], set(), None
    step_rows, long_rows, signals = [], [], []
    for k in range(a.max_signals + 1):
        res = {c: eng[c].run([m for m in cond if m in eng[c].order]) for c in cohorts}
        P = pd.DataFrame({c: res[c]['p'].reindex(universe) for c in cohorts})
        OR = pd.DataFrame({c: res[c]['or'].reindex(universe) for c in cohorts})
        if k == 0:
            marg = (P.copy(), OR.copy())
        sig = pd.DataFrame({c: P[c].lt(thr[c]['gec']) for c in cohorts})
        sig_any = set(sig.index[sig.any(axis=1)])
        if k > 0:
            lost = sorted(prev_sig - sig_any, key=lambda n: V.sort_key(n, a.cls))
            signals[-1]['members'] = ';'.join(lost)
            signals[-1]['n_members'] = len(lost)
        cand = [m for m in universe if m in sig_any and m not in cond]
        maxp = P.fillna(1.0).max(axis=1)
        pool = cand or [m for m in universe if m not in cond]
        lead = maxp.loc[pool].idxmin()
        entered = bool(cand) and k < a.max_signals
        reason = 'entered' if entered else ('step limit' if cand else 'no marker significant at GEC in any cohort')
        lp, lo = P.loc[lead].tolist(), OR.loc[lead].tolist()
        rec = dict(key, step=k, conditioned_on=';'.join(cond) or '-', lead=lead, entered=entered, reason=reason,
                   lead_max_p=float(maxp[lead]), n_sig_gec_any=len(sig_any - set(cond)))
        for s, c, x, o in zip(short, cohorts, lp, lo):
            rec[f'p_{s}'], rec[f'or_{s}'] = x, o
            rec[f'n_sig_gec_{s}'] = int(sig[c].sum())
        for t in V.THRESHOLDS:
            rec[t] = V.consistency(consistency(lp, lo, [thr[c][t] for c in cohorts]))
        step_rows.append(rec)
        for c in cohorts:
            r = res[c].reindex(members_of[c])           # an engine skips a marker with no minor allele
            r['errcode'] = r['errcode'].fillna('NO_RESULT')
            for m in members_of[c]:
                long_rows.append(dict(key, step=k, cohort=c, marker=m, p=r.at[m, 'p'], OR=r.at[m, 'or'],
                                      l95=r.at[m, 'l95'], u95=r.at[m, 'u95'], errcode=r.at[m, 'errcode'],
                                      sig_gec=bool(pd.notna(r.at[m, 'p']) and r.at[m, 'p'] < thr[c]['gec'])))
        if entered:
            mp, mo = marg[0].loc[lead].tolist(), marg[1].loc[lead].tolist()
            sg = dict(key, signal=len(signals) + 1, lead=lead, step=k, conditioned_on=';'.join(cond) or '-')
            for t in V.THRESHOLDS:
                sg[t] = rec[t]
            sg['gec_marginal'] = V.consistency(consistency(mp, mo, [thr[c]['gec'] for c in cohorts]))
            for s, x, o, xm, om in zip(short, lp, lo, mp, mo):
                sg[f'p_{s}'], sg[f'or_{s}'], sg[f'p_marginal_{s}'], sg[f'or_marginal_{s}'] = x, o, xm, om
            sg['carrier_alleles'] = ', '.join(carriers(lead)) if carriers else ''
            signals.append(sg)
            cond = cond + [lead]
        prev_sig = sig_any
        if not entered:
            break

    cls = a.cls
    pd.DataFrame(step_rows).to_csv(out / f'steps.{cls}.tsv', sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(long_rows).rename(columns={'OR': 'or'}).to_csv(out / f'stepwise.{cls}.tsv', sep='\t',
                                                                 index=False, float_format='%.6g')
    cols = list(key) + ['signal', 'lead', 'step', 'conditioned_on'] + list(V.THRESHOLDS) + ['gec_marginal'] + \
        [f'{w}_{s}' for s in short for w in ('p', 'or', 'p_marginal', 'or_marginal')] + \
        ['n_members', 'members', 'carrier_alleles']
    pd.DataFrame(signals, columns=cols).to_csv(out / f'signals.{cls}.tsv', sep='\t', index=False,
                                               float_format='%.6g')
    Path(out / f'leads.{cls}.txt').write_text(''.join(f'{s["lead"]}\n' for s in signals))
    runs = sum(e.n_runs for e in eng.values())
    print(f"[stepwise] {a.case_set} {a.model} {cls}: {len(signals)} signal(s) "
          f"{[(s['lead'], s['gec']) for s in signals]}; {len(step_rows)} steps, {runs} engine runs")


if __name__ == '__main__':
    main()
