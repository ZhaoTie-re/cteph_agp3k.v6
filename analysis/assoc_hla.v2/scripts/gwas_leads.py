#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The genome-wide random-model scan's MHC-region signals that this
#           component relates to the HLA markers, chosen by one stated rule.
#
#           Source: analysis/assoc_saige ('main' model: SAIGE, SEX + 10
#           bbj_mainland PCs, P = p.value.NA), its per-cohort peaks.tsv and
#           sumstats.tsv. Every peak whose span overlaps the MHC window is a
#           candidate; its lead variant is KEPT if its P is below the threshold
#           (1e-5, the scan's suggestive line) in EVERY cohort, read from each
#           cohort's own sumstats -- the three-cohort rule the project uses for
#           calling a signal. Every candidate is written with its per-cohort P
#           and the reason it was kept or excluded, so no MHC peak disappears
#           silently.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohorts', required=True, help='comma-separated')
    p.add_argument('--peaks', nargs='+', required=True, help="each cohort's peaks.tsv, in --cohorts order")
    p.add_argument('--sumstats', nargs='+', required=True, help="each cohort's sumstats.tsv, in --cohorts order")
    p.add_argument('--chrom', default='6')
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--end', type=int, required=True)
    p.add_argument('--p', type=float, default=1e-5)
    p.add_argument('--out', required=True)
    p.add_argument('--out-ids', required=True)
    return p.parse_args()


def lookup(sumstats, ids):
    """P, OR, A1 and A1 frequency of each ID in one sumstats file (grep, not a full read)."""
    pat = '\n'.join(f'\t{i}\t' for i in ids)
    r = subprocess.run(['grep', '-F', '-f', '/dev/stdin', sumstats], input=pat, capture_output=True, text=True)
    hdr = open(sumstats).readline().rstrip('\n').lstrip('#').split('\t')
    rows = [ln.split('\t') for ln in r.stdout.splitlines() if ln]
    d = pd.DataFrame(rows, columns=hdr) if rows else pd.DataFrame(columns=hdr)
    return d.set_index('ID')


def main():
    a = parse_args()
    cohorts = a.cohorts.split(',')
    if not (len(cohorts) == len(a.peaks) == len(a.sumstats)):
        sys.exit('ERROR: --cohorts, --peaks and --sumstats differ in length')
    cand = []
    for c, pk in zip(cohorts, a.peaks):
        p = pd.read_csv(pk, sep='\t', dtype={'chrom': str})
        hit = p[(p.chrom == a.chrom) & (p.end >= a.start) & (p.start <= a.end)]
        for _, r in hit.iterrows():
            cand.append({'lead_id': r['lead_id'], 'lead_pos': int(r['lead_pos']), 'found_in': c,
                         'peak_id': r['peak_id']})
    if not cand:
        sys.exit('ERROR: no genome-wide peak overlaps the MHC window; nothing to relate')
    cd = pd.DataFrame(cand)
    ids = sorted(set(cd.lead_id))
    rows = []
    for i in ids:
        sub = cd[cd.lead_id == i]
        rec = {'lead_id': i, 'pos': int(sub.lead_pos.iloc[0]),
               'peak_in': ','.join(f'{c}:{pid}' for c, pid in zip(sub.found_in, sub.peak_id))}
        ok_all = True
        for c, ss in zip(cohorts, a.sumstats):
            d = lookup(ss, [i])
            if i in d.index:
                rec[f'p_{c}'] = float(d.at[i, 'P'])
                rec[f'or_{c}'] = float(d.at[i, 'OR'])
                rec[f'a1_{c}'] = d.at[i, 'A1']
                rec[f'a1_freq_{c}'] = float(d.at[i, 'A1_FREQ'])
                ok_all &= rec[f'p_{c}'] < a.p
            else:
                rec[f'p_{c}'] = float('nan')
                ok_all = False
        rec['status'] = 'kept' if ok_all else 'excluded'
        failing = [c for c in cohorts if not rec.get(f'p_{c}', 1) < a.p]
        rec['reason'] = (f'P < {a.p:g} in all {len(cohorts)} cohorts' if ok_all
                         else f"P >= {a.p:g} (or absent) in {', '.join(failing)}")
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values('pos')
    out.to_csv(a.out, sep='\t', index=False, float_format='%.6g')
    kept = out[out.status == 'kept']
    if kept.empty:
        sys.exit(f'ERROR: no MHC lead has P < {a.p:g} in every cohort (see {a.out})')
    kept[['lead_id']].to_csv(a.out_ids, sep='\t', index=False, header=False)
    for _, r in out.iterrows():
        print(f"[gwas_leads] {r.lead_id}: {r.status} ({r.reason})")


if __name__ == '__main__':
    main()
