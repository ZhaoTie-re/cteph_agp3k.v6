#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Join every per-step table into the long table (one row per variant
#           x cohort) and run the checks that tie it to the upstream component:
#           recomputed lead statistics must equal lead_crosscohort.tsv, credible-
#           set rows must equal cs_variants_all.tsv, and the external frequencies
#           must agree with the study controls and with each other.
# Inputs  : outputs of select_signals / cs_status / genotype_stats /
#           annotate_function / ref_frequencies / opentargets; the upstream run
#           manifest; dbSNP (ANNOVAR avsnp table) to attribute rsIDs at
#           multi-allelic positions
# Outputs : saige_signals.long.tsv, checks.json
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import (as_bool, csv_list, load_json, load_variants, log, num,  # noqa: E402
                        read_table, read_tsv, write_tsv)

AF_COLS = ['af_cteph', 'af_agp3k', 'tommo_af', 'kg_jpt_af', 'kg_eas_af', 'gnomad_eas_af', 'gnomad_nfe_af']


OT_COLS = ['ot_in_index', 'ot_rsids', 'ot_clinvar_n', 'ot_clinvar', 'ot_gwas_n', 'ot_gwas_focus_n', 'ot_gwas_top',
           'ot_l2g', 'ot_qtl_n', 'ot_qtl_trans_n', 'ot_qtl_top', 'ot_e2g_n', 'ot_e2g_top']


# ── rsIDs at multi-allelic positions ─────────────────────────────────────────
# The upstream rsID comes from the ToMMo ID field, which at a multi-allelic position lists every rsID of the
# position, so it is withheld (`rsID_unresolved`). dbSNP is allele-specific: ANNOVAR's avsnp table, read
# through its index (100-bp bins with byte offsets), attributes the accession.
def annovar_allele(vid):
    """chr:pos:REF:ALT (VCF, anchored indels) -> (chrom, start, end, ref, alt) as ANNOVAR writes it."""
    c, p, ref, alt = vid.split(':')
    p, k = int(p), 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    r, a = ref[k:], alt[k:]
    if not r:                                         # insertion: after the last shared base
        return c.replace('chr', ''), p + k - 1, p + k - 1, '-', a
    return c.replace('chr', ''), p + k, p + k + len(r) - 1, r, a or '-'


def dbsnp_records(table, chrom, start, end):
    """dbSNP records (start, end, ref, alt, rsID) starting within [start, end] on chrom."""
    bins = sorted({b // 100 * 100 for b in range(start, end + 1)})
    pat = '^(' + '|'.join(f'{chrom}\t{b}' for b in bins) + ')\t'
    hits = subprocess.run(['grep', '-P', pat, f'{table}.idx'], capture_output=True, text=True).stdout
    out = []
    with open(table, 'rb') as fh:
        for line in hits.splitlines():
            _c, _b, s, e = line.split('\t')[:4]
            fh.seek(int(s))
            for rec in fh.read(int(e) - int(s)).decode().splitlines():
                f = rec.split('\t')
                if len(f) >= 6 and f[0] == chrom and start <= int(f[1]) <= end:
                    out.append((int(f[1]), int(f[2]), f[3], f[4], f[5]))
    return out


def resolve_rsids(v, table):
    """Attribute an rsID to a variant whose ToMMo ID field holds several accessions:
    1. the candidate dbSNP lists for this exact allele;
    2. else the only candidate that dbSNP does not list for another allele at the position.
    Otherwise the variant stays unresolved. Returns the updated table and an audit dict."""
    v = v.copy()
    v['rsID_source'] = np.where(v.rsID.fillna('.') != '.', 'ToMMo', '')
    audit = {}
    todo = v[(v.rsID.fillna('.') == '.') & (v.rsID_unresolved.fillna('.') != '.')]
    for i, r in todo.iterrows():
        cands = sorted({x for x in re.split(r'[;,/]', r.rsID_unresolved) if x.startswith('rs')})
        allele = annovar_allele(r.variant_id)
        chrom, start, end = allele[:3]
        recs = dbsnp_records(table, chrom, start - 1, end + 1) if table else []
        mine = {x[4] for x in recs if x[:4] == allele[1:] and x[4] in cands}
        other = {x[4] for x in recs if x[:4] != allele[1:] and x[4] in cands}
        pick, how = None, ''
        if len(mine) == 1:
            pick, how = mine.pop(), 'dbSNP allele'
        elif not mine and len(set(cands) - other) == 1:
            pick, how = (set(cands) - other).pop(), 'dbSNP exclusion'
        audit[r.variant_id] = {'candidates': cands, 'rsID': pick, 'rule': how or 'unresolved',
                               'dbsnp_records': [f'{s}-{e} {rf}>{al} {rs}' for s, e, rf, al, rs in recs
                                                 if rs in cands]}
        if pick:
            assert pick in cands, (r.variant_id, pick, cands)
            v.loc[i, ['rsID', 'rsID_unresolved', 'rsID_source']] = [pick, '.', how]
    return v, audit


def assemble_long(v, sel, status, geno, func, tommo, kg, gnomad, ot, cohorts, af_cohort):
    vv = v.merge(func, on='variant_id').merge(tommo, on='variant_id').merge(kg, on='variant_id') \
          .merge(gnomad, on='variant_id').merge(ot, on='variant_id')
    assert len(vv) == len(v)
    ref = geno[geno.cohort == af_cohort][['variant_id', 'case_eaf', 'ctrl_eaf']] \
        .rename(columns={'case_eaf': 'af_cteph', 'ctrl_eaf': 'af_agp3k'})
    vv = vv.merge(ref, on='variant_id')
    assert len(vv) == len(v)

    long = status.drop(columns=['locus_id']).merge(geno, on=['variant_id', 'cohort']).merge(vv, on='variant_id')
    assert len(long) == len(v) * len(cohorts)
    long['Beta'] = np.log(long.OR)
    # called_peak for the lead rows, from lead_crosscohort
    cp = {(r.variant_id, r.cohort): r.called_peak for _, r in sel.iterrows()}
    long['called_peak'] = [cp.get((a, b), '') for a, b in zip(long.variant_id, long.cohort)]

    notes = []
    for _, r in long.iterrows():
        n = []
        if pd.isna(r.tommo_af):
            n.append('ToMMo: not in release')
        elif r.tommo_filter != 'PASS':
            n.append(f'ToMMo: FILTER={r.tommo_filter}')
        if pd.isna(r.kg_eas_af):
            n.append('1KG: not in NYGC call set')
        if pd.isna(r.gnomad_eas_af) and pd.isna(r.gnomad_nfe_af):
            n.append('gnomAD: not present')
        notes.append('; '.join(n))
    long['external_af_note'] = notes
    rank = {c: i for i, c in enumerate(cohorts)}
    long['_c'] = long.cohort.map(rank)
    long = long.sort_values(['order', 'role_rank', 'max_pip', 'pos', '_c'],
                            ascending=[True, True, False, True, True]).drop(columns='_c')
    cols = ['order', 'locus_id', 'label', 'panel', 'tier', 'role', 'role_rank', 'variant_id', 'chrom', 'pos', 'EA', 'OA',
            'rsID', 'rsID_unresolved', 'rsID_source', 'Gene_EnsDb', 'Gene_Biotype', 'Gene_Distance_bp',
            'snpEff_Effect', 'snpEff_Impact', 'snpEff_Gene', 'snpEff_FeatureID', 'snpEff_Biotype',
            'snpEff_HGVS_c', 'snpEff_HGVS_p', 'snpEff_LoF_NMD', 'snpEff_n_annotations', 'snpEff_all_effects',
            'max_pip', 'cs_cohorts',
            'cohort', 'peak_id', 'peak_tier', 'called_peak', 'fine_mapped', 'fm_reason', 'cohort_has_cs', 'in_cs',
            'pip', 'pip_fitted', 'n_cs', 'cs_index', 'cs_size', 'cs_cut', 'cs_cuts_all', 'cs_status', 'cs_reason',
            'cs_reason_short', 'n_effects', 'n_sets_pre', 'n_sets_pure', 'dropped_set_size', 'dropped_set_min_r',
            'alpha_pure', 'alpha_cut_pure',
            'P', 'OR', 'L95', 'U95', 'Beta', 'SE', 'A1_FREQ', 'OBS_CT', 'n_case', 'n_ctrl',
            'case_homOA', 'case_het', 'case_homEA', 'case_missing_n', 'case_geno', 'case_eaf',
            'case_missing_rate', 'case_hwe_p',
            'ctrl_homOA', 'ctrl_het', 'ctrl_homEA', 'ctrl_missing_n', 'ctrl_geno', 'ctrl_eaf',
            'ctrl_missing_rate', 'ctrl_hwe_p',
            'af_cteph', 'af_agp3k', 'tommo_af', 'tommo_an', 'tommo_ac', 'tommo_filter', 'tommo_rsid',
            'tommo_match', 'tommo60_af', 'kg_eas_af', 'kg_eas_an', 'kg_jpt_af', 'kg_jpt_an', 'kg_filter',
            'kg_match', 'kg_eas_af_plink', 'kg_jpt_af_plink', 'gnomad_af', 'gnomad_eas_af', 'gnomad_nfe_af',
            'external_af_note'] + OT_COLS
    return long[cols].reset_index(drop=True)


def run_checks(long, sel, cs, min_tommo_frac):
    """The recomputed values must reproduce the component's own tables."""
    # leads: genotype string, EAF, missing rate, HWE P, OR / CI / P, group sizes
    m = long.merge(sel, on=['variant_id', 'cohort'], suffixes=('', '_ref'))
    assert len(m) == len(sel)
    for grp, pre in (('Case', 'case'), ('Control', 'ctrl')):
        geno_ref = m[f'{grp}_Genotype_Distribution']
        geno_new = m[f'{pre}_homOA'].astype(str) + '/' + m[f'{pre}_het'].astype(str) + '/' + m[f'{pre}_homEA'].astype(str)
        assert (geno_ref == geno_new).all(), m.loc[geno_ref != geno_new, ['variant_id', 'cohort']].to_string()
        assert (m[f'{pre}_eaf'].round(5) - num(m[f'{grp}_EAF'])).abs().max() < 1e-5
        assert (m[f'{pre}_missing_rate'].round(5) - num(m[f'{grp}_Missing_Rate'])).abs().max() < 1e-5
        assert (m[f'{pre}_hwe_p'] - num(m[f'{grp}_HWE_P'])).abs().max() < 1e-6
    assert np.allclose(m.P, num(m.P_ref), rtol=1e-6)
    for c in ('OR', 'L95', 'U95'):
        assert (m[c].round(5) - num(m[f'{c}_ref'])).abs().max() < 1e-5, c
    assert (m.n_case == num(m.N_case)).all() and (m.n_ctrl == num(m.N_ctrl)).all()
    # credible-set rows: P / OR / PIP (cs_variants_all rounds PIP to 6 dp)
    mc = long.merge(cs, on=['variant_id', 'cohort'], suffixes=('', '_ref'))
    assert len(mc) == len(cs) and mc.in_cs.all()
    assert np.allclose(mc.P, num(mc.P_ref), rtol=1e-6) and np.allclose(mc.pip, num(mc.pip_num), atol=1e-6)
    assert np.allclose(mc.OR, num(mc.OR_ref), atol=1e-5)
    # external frequencies
    vals = long.drop_duplicates('variant_id')[AF_COLS]
    assert ((vals.fillna(0) >= 0) & (vals.fillna(0) <= 1)).all().all()
    assert vals.notna().any(axis=1).all()
    both = vals.dropna(subset=['tommo_af', 'af_agp3k'])
    r1 = np.corrcoef(both.tommo_af, both.af_agp3k)[0, 1]
    both2 = vals.dropna(subset=['gnomad_eas_af', 'kg_eas_af'])
    r2 = np.corrcoef(both2.gnomad_eas_af, both2.kg_eas_af)[0, 1]
    log(f'r(ToMMo, AGP3K EAF) = {r1:.4f} on {len(both)}; r(gnomAD EAS, 1KG EAS) = {r2:.4f} on {len(both2)}')
    assert r1 > 0.98 and r2 > 0.95
    n_v = len(vals)
    assert vals.tommo_af.notna().sum() >= min_tommo_frac * n_v
    miss_any = long[AF_COLS[2:]].isna().any(axis=1)
    assert (long.loc[miss_any, 'external_af_note'] != '').all()
    v1 = long.drop_duplicates('variant_id')
    return {
        'n_variants': n_v, 'n_rows': len(long),
        'n_lead_rows_checked': len(m), 'n_cs_rows_checked': len(mc),
        'r_tommo_agp3k': round(float(r1), 4), 'r_gnomad_eas_1kg_eas': round(float(r2), 4),
        'coverage': {'snpeff': int((v1.snpEff_Effect != '.').sum()), 'tommo': int(v1.tommo_af.notna().sum()),
                     'kg_eas': int(v1.kg_eas_af.notna().sum()), 'kg_jpt': int(v1.kg_jpt_af.notna().sum()),
                     'gnomad_eas': int(v1.gnomad_eas_af.notna().sum()),
                     'gnomad_nfe': int(v1.gnomad_nfe_af.notna().sum())},
        'tommo_absent': sorted(v1[v1.tommo_match == 'absent'].variant_id),
        'tommo_normalised': sorted(v1[v1.tommo_match == 'normalised'].variant_id),
        'tommo_filtered': sorted(f'{r.variant_id} ({r.tommo_filter})' for _, r in
                                 v1[(v1.tommo_match != 'absent') & (v1.tommo_filter != 'PASS')].iterrows()),
        'kg_absent': sorted(v1[v1.kg_match == 'absent'].variant_id),
        'lof_nmd': int((v1.snpEff_LoF_NMD.fillna('') != '').sum()),
        'tommo_an_median': int(v1.tommo_an.dropna().median()) if v1.tommo_an.notna().any() else None,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--loci', required=True)
    ap.add_argument('--sel', required=True)
    ap.add_argument('--cs-members', required=True)
    ap.add_argument('--selection-json', required=True)
    ap.add_argument('--status', required=True)
    ap.add_argument('--cs-json', required=True)
    ap.add_argument('--geno', nargs='+', required=True)
    ap.add_argument('--function', required=True)
    ap.add_argument('--tommo', required=True)
    ap.add_argument('--kg', required=True)
    ap.add_argument('--gnomad', required=True)
    ap.add_argument('--opentargets', required=True, help='opentargets.tsv (one row per variant)')
    ap.add_argument('--dbsnp', default='', help='ANNOVAR avsnp table (with .idx) to attribute rsIDs; empty: skip')
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--af-cohort', required=True)
    ap.add_argument('--upstream-manifest', required=True)
    ap.add_argument('--p-value', required=True, choices=['spa', 'no_spa'])
    ap.add_argument('--min-tommo-frac', type=float, default=0.9)
    ap.add_argument('--out-long', required=True)
    ap.add_argument('--out-checks', required=True)
    a = ap.parse_args()
    cohorts = csv_list(a.cohorts)
    assert a.af_cohort in cohorts
    up = load_json(a.upstream_manifest)
    assert up['p_value_source'] == a.p_value, (up['p_value_source'], a.p_value)

    v = load_variants(a.variants)
    v, rsid_audit = resolve_rsids(v, a.dbsnp)
    for vid, x in rsid_audit.items():
        log(f'rsID {vid}: {"/".join(x["candidates"])} -> {x["rsID"] or "unresolved"} ({x["rule"]}; '
            f'{"; ".join(x["dbsnp_records"]) or "no dbSNP record of the candidates"})')
    sel = read_tsv(a.sel)
    cs = read_tsv(a.cs_members)
    status = read_table(a.status)
    for c in ('peak_id', 'peak_tier', 'fm_reason', 'cs_cuts_all', 'cs_reason', 'cs_reason_short'):
        status[c] = status[c].fillna('').astype(str)
    for c in ('in_cs', 'fine_mapped', 'pip_fitted', 'cohort_has_cs'):
        status[c] = as_bool(status[c])
    geno = pd.concat([read_table(p) for p in a.geno], ignore_index=True)
    assert set(geno.cohort) == set(cohorts) and len(geno) == len(v) * len(cohorts)

    def ref(path, text):
        d = read_table(path)
        for c in text:
            if c in d.columns:
                d[c] = d[c].fillna('').astype(str)
        return d
    func = ref(a.function, ['snpEff_Effect', 'snpEff_Impact', 'snpEff_Gene', 'snpEff_FeatureID', 'snpEff_Biotype',
                            'snpEff_HGVS_c', 'snpEff_HGVS_p', 'snpEff_all_effects', 'snpEff_LoF_NMD',
                            'Gene_EnsDb', 'Gene_Biotype'])
    tommo = ref(a.tommo, ['tommo_filter', 'tommo_rsid', 'tommo_match'])
    kg = ref(a.kg, ['kg_filter', 'kg_match'])
    gnomad = read_table(a.gnomad)
    ot = ref(a.opentargets, ['ot_rsids', 'ot_clinvar', 'ot_gwas_top', 'ot_l2g', 'ot_qtl_top', 'ot_e2g_top'])
    for c in OT_COLS:
        if c not in ot.columns:                      # Open Targets not queried: empty columns
            ot[c] = '' if c in ('ot_rsids', 'ot_clinvar', 'ot_gwas_top', 'ot_l2g', 'ot_qtl_top', 'ot_e2g_top') else np.nan
    assert ot.variant_id.is_unique and set(ot.variant_id) == set(v.variant_id)
    counts = ['ot_clinvar_n', 'ot_gwas_n', 'ot_gwas_focus_n', 'ot_qtl_n', 'ot_qtl_trans_n', 'ot_e2g_n']
    if ot.ot_in_index.notna().all():
        ot['ot_in_index'] = as_bool(ot.ot_in_index)
        assert (ot[counts] >= 0).all().all()
        assert (ot.ot_qtl_trans_n <= ot.ot_qtl_n).all() and (ot.ot_gwas_focus_n <= ot.ot_gwas_n).all()
        assert (ot.loc[~ot.ot_in_index, counts] == 0).all().all(), 'a variant absent from the index has records'

    long = assemble_long(v, sel, status, geno, func, tommo, kg, gnomad, ot[['variant_id'] + OT_COLS], cohorts,
                         a.af_cohort)
    checks = run_checks(long, sel, cs, a.min_tommo_frac)
    v1 = long.drop_duplicates('variant_id')
    if v1.ot_in_index.notna().all():
        checks['opentargets'] = {'in_index': int(v1.ot_in_index.astype(bool).sum()),
                                 'with_gwas': int((v1.ot_gwas_n > 0).sum()),
                                 'with_thrombosis_gwas': int((v1.ot_gwas_focus_n > 0).sum()),
                                 'with_molqtl': int((v1.ot_qtl_n > 0).sum()),
                                 'with_trans_molqtl': int((v1.ot_qtl_trans_n > 0).sum()),
                                 'with_e2g': int((v1.ot_e2g_n > 0).sum()),
                                 'with_clinvar': int((v1.ot_clinvar_n > 0).sum())}
    checks['rsid_resolution'] = {'dbsnp': Path(a.dbsnp).name if a.dbsnp else '', 'variants': rsid_audit,
                                 'unresolved': sorted(v1[v1.rsID_unresolved.fillna('.') != '.'].variant_id)}
    selection = load_json(a.selection_json)
    cs_json = load_json(a.cs_json)
    checks = {'all_passed': True, 'built': dt.datetime.now().isoformat(timespec='seconds'),
              'p_value_source': up['p_value_source'], 'upstream_p_value_column': up['p_value_column'],
              'cohorts': cohorts, 'af_cohort': a.af_cohort,
              'selection': selection, 'credible_sets': cs_json, 'table': checks}
    write_tsv(long, a.out_long)
    with open(a.out_checks, 'w') as fh:
        json.dump(checks, fh, indent=2)
    log(f'long table: {len(long)} rows; CHECKS PASSED')


if __name__ == '__main__':
    main()
