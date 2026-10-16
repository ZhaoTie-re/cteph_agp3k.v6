#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Per-cohort association statistics and SuSiE PIP for every table
#           variant x cohort. Inside a cohort's fine-mapping window the values
#           come from that peak's own pip.tsv / sumstat.tsv; a variant outside
#           the window, or a locus with no peak in the cohort, falls back to the
#           cohort's genome-wide sumstats (no PIP: "not fine-mapped"). The
#           genome-wide file is read for every variant so the two sources are
#           asserted to agree wherever both exist.
# Inputs  : variants.tsv, pairs.tsv, cs_members.tsv (select_signals.py);
#           --peak KEY PIP SUMSTAT per locus x cohort with a peak (KEY = cohort/peak_id);
#           --gwas cohort=<cohort>/01.assoc/<peak model>/sumstats.tsv
# Outputs : assoc.tsv  one row per variant x cohort
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import (as_bool, csv_list, kv_paths, load_variants, log, num,  # noqa: E402
                        read_tsv, run, write_tsv)


def gwas_rows(path, ids, tag):
    """Rows of a genome-wide sumstats.tsv for a set of IDs (awk hash filter; files are ~0.7 GB)."""
    head = open(path).readline().rstrip('\n').split('\t')
    assert head[2] == 'ID', head
    idfile = Path(f'ids.{tag}.txt')
    idfile.write_text('\n'.join(sorted(ids)) + '\n')
    out = run(['awk', '-F\t', 'NR==FNR{w[$1];next} FNR==1||($3 in w)', idfile, path])
    idfile.unlink()
    rows = [l.split('\t') for l in out.splitlines()]
    g = pd.DataFrame(rows[1:], columns=rows[0]).rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
    assert (g.A1 == g.ALT).all(), 'A1 is not ALT in the genome-wide sumstats'
    return g.set_index('ID')


def fetch_assoc(v, pairs, cs, peak_files, gwas):
    gw = {c: gwas_rows(p, set(v.variant_id), c) for c, p in gwas.items()}
    rows = []
    for _, pr in pairs.iterrows():
        locus, cohort, peak, tier = pr.locus_id, pr.cohort, pr.peak_id, pr.peak_tier
        want = v[v.locus_id == locus]
        m = pd.DataFrame(columns=['ID'])
        set_ids, sizes, cuts, pip_fitted = [], {}, {}, False
        if peak:
            pip_path, ss_path = peak_files[f'{cohort}/{peak}']
            pip = read_tsv(pip_path)
            ss = read_tsv(ss_path)
            pip_num = num(pip.pip)
            pip_fitted = bool((pip_num > 0).any())          # SuSiE can return an all-zero PIP vector
            in_set = ~pip.cs.isin(['NA', '', '.'])
            set_ids = sorted(int(float(x)) for x in pip.cs[in_set].unique())      # 1..k, k = 0 without a set
            sizes = {k: int((pip.cs == str(k)).sum()) for k in set_ids}
            cuts = {k: float(pip_num[pip.cs == str(k)].min()) for k in set_ids}   # smallest PIP inside set k
            m = pip.merge(ss[['ID', 'REF', 'ALT', 'L95', 'U95', 'OBS_CT']], on='ID', how='left')
            m = m[m.ID.isin(want.variant_id)]
            assert (m.A1 == m.ALT).all(), (locus, cohort)
        has_cs = (locus, cohort) in {(l, c) for l, c in zip(cs.locus_id, cs.cohort)}
        assert has_cs == (len(set_ids) > 0), (locus, cohort)
        common = {'cohort': cohort, 'peak_id': peak, 'peak_tier': tier, 'cohort_has_cs': has_cs,
                  'n_cs': len(set_ids), 'cs_cuts_all': ';'.join(f'{k}:{cuts[k]:.6g}' for k in set_ids)}
        for _, r in m.iterrows():
            in_cs = r.cs not in ('NA', '', '.')
            k = int(float(r.cs)) if in_cs else np.nan
            rows.append({**common, 'variant_id': r.ID, 'fine_mapped': True, 'fm_reason': '',
                         'P': float(r.P), 'OR': float(r.OR), 'SE': float(r.SE),
                         'L95': float(r.L95), 'U95': float(r.U95),
                         'A1_FREQ': float(r.A1_FREQ), 'OBS_CT': int(float(r.OBS_CT)),
                         'pip': float(r.pip), 'pip_fitted': pip_fitted, 'in_cs': in_cs, 'cs_index': k,
                         'cs_size': sizes.get(k, np.nan) if in_cs else np.nan,
                         'cs_cut': cuts.get(k, np.nan) if in_cs else (min(cuts.values()) if cuts else np.nan)})
        # the fallback: not in this cohort's fine-mapping window, or no peak here
        for vid in sorted(set(want.variant_id) - set(m.ID)):
            reason = 'outside_window' if peak else 'no_peak'
            rec = {**common, 'variant_id': vid, 'fine_mapped': False, 'pip': np.nan, 'pip_fitted': False,
                   'in_cs': False, 'cs_index': np.nan, 'cs_size': np.nan, 'cs_cut': np.nan}
            if vid in gw[cohort].index:
                g = gw[cohort].loc[vid]
                rec.update({'fm_reason': reason, 'P': float(g.P), 'OR': float(g.OR), 'SE': float(g.SE),
                            'L95': float(g.L95), 'U95': float(g.U95), 'A1_FREQ': float(g.A1_FREQ),
                            'OBS_CT': int(float(g.OBS_CT))})
            else:
                rec.update({'fm_reason': 'not_tested', 'P': np.nan, 'OR': np.nan, 'SE': np.nan,
                            'L95': np.nan, 'U95': np.nan, 'A1_FREQ': np.nan, 'OBS_CT': np.nan})
            rows.append(rec)
    a = pd.DataFrame(rows)
    assert len(a) == len(v) * pairs.cohort.nunique(), (len(a), len(v))
    assert not a.duplicated(['variant_id', 'cohort']).any()
    tested = a.fm_reason != 'not_tested'
    assert a.loc[tested, ['P', 'OR', 'L95', 'U95']].notna().all().all()
    assert (a.pip.notna() == a.fine_mapped).all()
    # in_cs from pip.tsv must be the membership cs_variants_all reports
    mem = set(zip(cs.variant_id, cs.cohort))
    assert set(zip(a[a.in_cs].variant_id, a[a.in_cs].cohort)) == mem
    # the peak files and the genome-wide scan are the same numbers
    fm = a[a.fine_mapped]
    ref = np.array([float(gw[c].loc[vid].P) for vid, c in zip(fm.variant_id, fm.cohort)])
    assert np.allclose(fm.P.to_numpy(), ref, rtol=1e-6), 'peak-file P differs from the genome-wide sumstats'
    n_fb = int((~a.fine_mapped).sum())
    log(f'association: {len(a)} variant x cohort rows; {n_fb} not fine-mapped '
        f'({a[~a.fine_mapped].fm_reason.value_counts().to_dict() if n_fb else {}})')
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--cs-members', required=True)
    ap.add_argument('--peak', nargs=3, action='append', default=[], metavar=('KEY', 'PIP', 'SUMSTAT'))
    ap.add_argument('--gwas', action='append', required=True, help='cohort=<peak model>/sumstats.tsv')
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    cohorts = csv_list(a.cohorts)
    v = load_variants(a.variants)
    pairs = read_tsv(a.pairs)
    cs = read_tsv(a.cs_members)
    gwas = kv_paths(a.gwas)
    assert set(gwas) == set(cohorts)
    peak_files = {k: (Path(p), Path(s)) for k, p, s in a.peak}
    need = {f'{r.cohort}/{r.peak_id}' for _, r in pairs[pairs.peak_id != ''].iterrows()}
    assert need == set(peak_files), (need ^ set(peak_files))
    out = fetch_assoc(v, pairs, cs, peak_files, gwas)
    out['in_cs'] = as_bool(out.in_cs)
    write_tsv(out, a.out)


if __name__ == '__main__':
    main()
