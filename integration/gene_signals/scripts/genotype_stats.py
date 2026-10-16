#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One cohort: genotypes of every table variant on the call set the gene
#           tests ran on, and the carriers of every table gene.
#             variants  REF/REF, REF/ALT, ALT/ALT, missing per cases and controls;
#                       ALT frequency; missing rate; HWE P (exact, within the group);
#                       pooled minor allele and MAC;
#                       whether the cohort's test of the gene included it, and if
#                       not, why (MAC below the cohort's minAC; control MAF >= the
#                       rare-variant cut in this cohort, with that MAF; removed by the 15x-30x
#                       missingness-gap filter; not in the call set)
#             carriers  per gene x variant set: carriers of >= 1 minor allele among
#                       the variants the test used, by case platform and control
#                       source, multi-site carriers, cumulative MAC
#           The minor allele is decided once on the pooled counts, as the gene
#           tests do. Every count of a tested variant and every gene total is
#           asserted equal to the upstream tables.
# Inputs  : variants.tsv; the cohort's call set (PLINK); the variant-QC exclude
#           list; the common-variant list and the control frequencies behind it;
#           the sample sheet (platform, source);
#           upstream variant_stats_all.tsv and set_stats_all.tsv
# Outputs : geno.<cohort>.tsv (one row per variant), carriers.<cohort>.tsv (long), samples.<cohort>.tsv
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import csv_list, import_annot, load_variants, log, num, read_tsv, run, write_tsv  # noqa: E402


def in_list(path, ids, work):
    """The ids found in a one-ID-per-line list (grep -Fx: the lists hold millions of lines)."""
    f = work / f'ids.{Path(path).name}.txt'
    f.write_text('\n'.join(sorted(ids)) + '\n')
    try:
        out = run(['grep', '-Fx', '-f', f, path])
    except RuntimeError:                       # grep exits 1 when nothing matches
        return set()
    return set(out.split())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--cohort-code', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--bfile', required=True, help='PLINK prefix of the call set the gene tests ran on')
    ap.add_argument('--exclude-list', required=True, help='variants removed by the 15x-30x missingness-gap filter')
    ap.add_argument('--common-list', required=True, help='variants with control MAF >= the rare-variant cut')
    ap.add_argument('--control-afreq', required=True, help='plink2 --freq of the controls the rare/common split used')
    ap.add_argument('--rare-maf', type=float, required=True, help='the rare-variant cut on control MAF')
    ap.add_argument('--minac', type=int, required=True, help="the cohort's minor-allele-count floor")
    ap.add_argument('--sample-info', required=True, help='sample sheet: ID, WGS_Platform, Project')
    ap.add_argument('--variant-stats', required=True, help='upstream variant_stats_all.tsv (cross-check)')
    ap.add_argument('--set-stats', required=True, help='upstream set_stats_all.tsv (cross-check)')
    ap.add_argument('--plink2', required=True)
    ap.add_argument('--annot-lib', required=True, help='the staged _shared/scripts/variant_annot.py (HWE)')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--out-geno', required=True)
    ap.add_argument('--out-carriers', required=True)
    ap.add_argument('--out-samples', required=True, help='the cohort\'s samples: IID, group, platform, source')
    ap.add_argument('--out-multisite', required=True, help='multi-site carriers: one row per gene x set x sample '
                    '(sample IDs: analyst file)')
    a = ap.parse_args()
    cohort, code, strata = a.cohort, a.cohort_code, csv_list(a.strata)
    work = Path('work_geno')
    work.mkdir(exist_ok=True)
    v = load_variants(a.variants)
    ids = sorted(set(v.variant_id))

    # ── samples ────────────────────────────────────────────────────────────
    fam = pd.read_csv(f'{a.bfile}.fam', sep=r'\s+', header=None, names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'],
                      dtype={'FID': str, 'IID': str})
    assert set(fam.PHENO) == {1, 2}, set(fam.PHENO)
    info = pd.read_csv(a.sample_info, dtype=str).set_index('ID_JHRPv6')
    meta = info.reindex(fam.IID)
    assert meta.WGS_Platform.notna().all() and meta.Project.notna().all(), 'samples missing from the sample sheet'
    is_case = (fam.PHENO == 2).to_numpy()
    n_case, n_ctrl = int(is_case.sum()), int((~is_case).sum())

    # ── genotypes of the variants in the call set ──────────────────────────
    bim = pd.read_csv(f'{a.bfile}.bim', sep='\t', header=None, usecols=[1], dtype=str)[1]
    present = sorted(set(ids) & set(bim))
    (work / 'extract.txt').write_text('\n'.join(present) + '\n')
    run([a.plink2, '--bfile', a.bfile, '--extract', work / 'extract.txt', '--export', 'A',
         '--threads', a.threads, '--out', work / 'g'])
    raw = pd.read_csv(work / 'g.raw', sep=r'\s+', dtype={'FID': str, 'IID': str})
    assert raw.IID.tolist() == fam.IID.tolist(), 'sample order differs from the .fam'
    alt = {}                                   # variant -> ALT count per sample (NaN = missing)
    for col in raw.columns[6:]:
        vid, counted = col.rsplit('_', 1)
        _c, _p, ref, al = vid.split(':')
        x = raw[col].to_numpy(dtype=float)
        assert counted in (ref, al), col
        alt[vid] = x if counted == al else 2 - x
    assert set(alt) == set(present), set(present) ^ set(alt)

    # ── why a variant is absent ────────────────────────────────────────────
    absent = sorted(set(ids) - set(present))
    excl = in_list(a.exclude_list, absent, work) if absent else set()
    common = in_list(a.common_list, absent, work) if absent else set()
    # the control MAF that put a variant on the common side, from the split's own frequency table
    common_maf = {}
    if common:
        (work / 'common.ids').write_text('\n'.join(sorted(common)) + '\n')
        out = run(['awk', '-F', '\t', 'NR==FNR{w[$1];next} ($2 in w){print $2"\t"$6}', work / 'common.ids',
                   a.control_afreq])
        for line in out.splitlines():
            vid, f = line.split('\t')
            common_maf[vid] = min(float(f), 1 - float(f))
        assert set(common_maf) == common, common ^ set(common_maf)
        assert all(m >= a.rare_maf for m in common_maf.values()), common_maf

    rows = []
    for vid in ids:
        rec = {'variant_id': vid, 'cohort': cohort, 'n_case': n_case, 'n_ctrl': n_ctrl}
        if vid not in alt:
            rec['in_callset'] = False
            rec['callset_reason'] = ('gap_filter' if vid in excl else 'common' if vid in common else 'not_in_callset')
            rec['ctrl_maf_split'] = common_maf.get(vid, np.nan)
            rows.append(rec)
            continue
        x = alt[vid]
        called = ~np.isnan(x)
        n_alt = np.nansum(x)
        n_obs = 2 * called.sum()
        minor_is_alt = n_alt <= n_obs - n_alt         # pooled, as the gene tests decide it
        m = np.where(called, x if minor_is_alt else 2 - x, 0)
        rec.update({'in_callset': True, 'callset_reason': '', 'minor_allele': 'ALT' if minor_is_alt else 'REF',
                    'mac': int(m.sum()), 'carrier_iids': ';'.join(fam.IID[m > 0])})
        for grp, sel in (('case', is_case), ('ctrl', ~is_case)):
            xs, cs = x[sel], called[sel]
            rec.update({f'{grp}_homRef': int(((xs == 0) & cs).sum()), f'{grp}_het': int(((xs == 1) & cs).sum()),
                        f'{grp}_homAlt': int(((xs == 2) & cs).sum()), f'{grp}_missing_n': int((~cs).sum()),
                        f'{grp}_called': int(cs.sum()),
                        f'{grp}_alt_af': float(np.nansum(xs) / (2 * cs.sum())) if cs.sum() else np.nan,
                        f'{grp}_missing_rate': float((~cs).sum() / sel.sum()),
                        f'{grp}_mac': int(m[sel].sum()), f'{grp}_carriers': int((m[sel] > 0).sum())})
        rows.append(rec)
    g = pd.DataFrame(rows)

    # ── HWE P per group: the shared group_stats (plink2 --geno-counts --hardy), the function behind the
    # common-variant table; its genotype counts must be the ones exported above
    va = import_annot(a.annot_lib)
    gi0 = g.set_index('variant_id')
    for grp, code_, n_grp in (('case', 2, n_case), ('ctrl', 1, n_ctrl)):
        st, n = va.group_stats(a.plink2, a.bfile, fam, work / 'hwe', grp, code_, present, threads=a.threads)
        assert n == n_grp and set(st) == set(present), (cohort, grp, n, len(st))   # {} when plink2 fails
        for vid, d in st.items():
            x = gi0.loc[vid]
            assert (d['ref'], d['alt']) == tuple(vid.split(':')[2:4]), (vid, d['ref'], d['alt'])
            assert (int(d['hom_ref']), int(d['het']), int(d['hom_alt'])) == \
                (int(x[f'{grp}_homRef']), int(x[f'{grp}_het']), int(x[f'{grp}_homAlt'])), (cohort, grp, vid)
        g[f'{grp}_hwe_p'] = g.variant_id.map({vid: float(d['hwe_p']) for vid, d in st.items()})

    # ── test membership: the variant sets that tested the variant in this cohort ──
    vv = v.drop_duplicates('variant_id').set_index('variant_id')
    for s in strata:
        tested = {vid for vid in ids if code in str(vv.loc[vid, f'cohorts_{s}'])}
        g[f'in_set_{s}'] = g.variant_id.isin(tested)
        reason = []
        for _, r in g.iterrows():
            if r[f'in_set_{s}']:
                reason.append('')
            elif not r.in_callset:
                reason.append(r.callset_reason)
            elif s == 'moderate_high' and vv.loc[r.variant_id, 'impact'] == 'LOW':
                reason.append('impact')
            elif r.mac < a.minac:
                reason.append('mac')
            else:
                reason.append('unexplained')
        g[f'set_reason_{s}'] = reason
        bad = g[g[f'set_reason_{s}'] == 'unexplained']
        assert bad.empty, f'{cohort}/{s}: variants in the call set with MAC >= {a.minac} but not tested: ' \
                          f'{bad.variant_id.tolist()}'
        bad = g[g[f'in_set_{s}'] & ((~g.in_callset) | (g.mac < a.minac))]
        assert bad.empty, f'{cohort}/{s}: tested variants below minAC or absent: {bad.variant_id.tolist()}'

    # ── cross-check every tested variant against the upstream counts ───────
    up = read_tsv(a.variant_stats)
    up = up[(up.cohort == cohort) & up.variant_id.isin(ids)]
    gi = g.set_index('variant_id')
    n_chk = 0
    for _, r in up.iterrows():
        x = gi.loc[r.variant_id]
        flipped = r.flipped == '1'
        assert (x.minor_allele == 'REF') == flipped, (r.variant_id, x.minor_allele, r.flipped)
        exp = {'case_mac': int(r.mac_case), 'ctrl_mac': int(r.mac_control),
               'case_carriers': int(r.n_carrier_case), 'ctrl_carriers': int(r.n_carrier_control),
               'case_missing_n': int(r.n_missing_case), 'ctrl_missing_n': int(r.n_missing_control)}
        got = {k: int(x[k]) for k in exp}
        assert got == exp, (cohort, r.variant_id, got, exp)
        n_chk += 1

    # ── carriers per gene x variant set ────────────────────────────────────
    plat, proj = meta.WGS_Platform.to_numpy(), meta.Project.to_numpy()
    ss = read_tsv(a.set_stats)
    ss = ss[ss.cohort == cohort].set_index(['set_name', 'stratum'])
    crows, mrows = [], []
    for sn, gv in v.groupby('set_name', sort=False):
        for s in strata:
            tested = [vid for vid in gv.variant_id if bool(gi.loc[vid, f'in_set_{s}'])]
            if not tested:
                continue
            M = np.vstack([np.where(np.isnan(alt[vid]), 0,
                                    alt[vid] if gi.loc[vid, 'minor_allele'] == 'ALT' else 2 - alt[vid])
                           for vid in tested])
            carrier = (M > 0).any(0)
            sites = (M > 0).sum(0)
            mac = M.sum(0)
            base = {'set_name': sn, 'stratum': s, 'cohort': cohort, 'n_var_tested': len(tested)}
            for grp, sel in (('case', is_case), ('ctrl', ~is_case)):
                crows.append(dict(base, group=grp, by='all', subgroup='all', n_samples=int(sel.sum()),
                                  carriers=int(carrier[sel].sum()), multi_site=int((sites[sel] >= 2).sum()),
                                  mac=int(mac[sel].sum())))
                lab = plat if grp == 'case' else proj
                for sub in sorted(set(lab[sel])):
                    k = sel & (lab == sub)
                    crows.append(dict(base, group=grp, by='platform' if grp == 'case' else 'source', subgroup=sub,
                                      n_samples=int(k.sum()), carriers=int(carrier[k].sum()),
                                      multi_site=int((sites[k] >= 2).sum()), mac=int(mac[k].sum())))
            # the totals must be the upstream set statistics
            u = ss.loc[(sn, s)]
            got = (int(carrier[is_case].sum()), int(carrier[~is_case].sum()), int(mac[is_case].sum()),
                   int(mac[~is_case].sum()), int((sites[is_case] >= 2).sum()), int((sites[~is_case] >= 2).sum()),
                   len(tested))
            exp = (int(u.n_carrier_case), int(u.n_carrier_control), int(u.mac_case), int(u.mac_control),
                   int(u.n_multi_case), int(u.n_multi_control), int(u.n_var))
            assert got == exp, (cohort, sn, s, got, exp)
            # the multi-site carriers themselves: the tested variants each carries a minor allele at
            for k in np.where(sites >= 2)[0]:
                carried = [(vid, int(M[q, k])) for q, vid in enumerate(tested) if M[q, k] > 0]
                mrows.append({'set_name': sn, 'stratum': s, 'cohort': cohort, 'IID': fam.IID.iloc[k],
                              'group': 'case' if is_case[k] else 'ctrl',
                              'platform': plat[k], 'source': proj[k], 'n_variants_carried': len(carried),
                              'n_hom': sum(d == 2 for _v, d in carried),
                              'variant_ids': ';'.join(vid for vid, _d in carried),
                              'dosages': ';'.join(str(d) for _v, d in carried)})
    car = pd.DataFrame(crows)
    ms = pd.DataFrame(mrows, columns=['set_name', 'stratum', 'cohort', 'IID', 'group', 'platform', 'source',
                                      'n_variants_carried', 'n_hom', 'variant_ids', 'dosages'])
    # one multi-site row per counted multi-site carrier
    tot_ms = car[car.by == 'all'].groupby(['set_name', 'stratum', 'group']).multi_site.sum()
    got_ms = ms.groupby(['set_name', 'stratum', 'group']).size()
    assert tot_ms[tot_ms > 0].sort_index().to_dict() == got_ms.sort_index().to_dict(), (cohort, 'multi-site rows')
    write_tsv(g, a.out_geno)
    pd.DataFrame({'IID': fam.IID, 'group': np.where(is_case, 'case', 'ctrl'), 'platform': plat,
                  'source': proj}).to_csv(a.out_samples, sep='\t', index=False)
    write_tsv(car, a.out_carriers)
    write_tsv(ms, a.out_multisite)
    log(f'{cohort}: {n_case} cases / {n_ctrl} controls; {len(present)}/{len(ids)} variants in the call set '
        f'(absent: {len(excl)} gap filter, {len(common)} common, {len(absent) - len(excl) - len(common)} other); '
        f'{n_chk} tested variants and {car[(car.by == "all") & (car.group == "case")].shape[0]} gene x set totals '
        f'equal the upstream counts')


if __name__ == '__main__':
    main()
