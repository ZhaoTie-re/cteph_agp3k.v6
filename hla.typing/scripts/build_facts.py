#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ONE SOURCE FOR EVERY NUMBER THAT APPEARS IN PROSE.
#
#           This component's figures, sidecars and documents used to quote
#           hand-typed values. When the frequency reference changed from the
#           1000 Genomes JPT panel to jMorp 61KJPN, the tables followed and the
#           prose did not. Captions still said "105 samples" and "Five loci"  [retired]
#           over a figure drawn from a far larger panel at 13 loci, and the caveat
#           "105 samples cannot sample a rare allele" was left arguing the       [retired]
#           opposite of what the data now support.
#
#           A literal cannot follow its data. So no number is written into prose
#           anywhere in this component; every one of them is derived here, from
#           the published tables, and rendered ONCE in the form prose uses.
#           Figures interpolate these values, documents quote them, the report
#           reads this file, and verify.sh fails on any number in a .md that is
#           not in here or in docs/NUMBERS_ALLOWED.txt.
#
#           Read-only over results/. It never touches the association's inputs.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLE_ID = {'sample_id': str}

# The level above which a locus is read at face value. It is a threshold we
# chose, not a result, but the prose quotes it, so it is published as a fact
# rather than typed into two report editions that could then drift apart.
R_CONCORDANT = 0.98


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', required=True, help='the component results/ directory')
    p.add_argument('--workbook', default=None,
                   help='the v6 sample workbook; the cohort denominator comes from it')
    p.add_argument('--id-col', default='ID_JHRPv6')
    p.add_argument('--flag-col', default='Flag_JHRPv6')
    p.add_argument('--group-col', default='Outcome')
    p.add_argument('--platform-col', default='WGS_Platform')
    p.add_argument('--keep-id', default=None, help="wgs.auto.par's sample_qc.keep.id")
    p.add_argument('--cohort-keep', default=None,
                   help='the analysis cohort list, e.g. full_mainland.fid_iid.txt')
    p.add_argument('--cohort-name', default='full_mainland')
    p.add_argument('--cohort-pheno', nargs='*', default=[],
                   help='NAME=path/to/pheno.tsv for each downstream cohort, to record containment')
    p.add_argument('--reference-jmorp', default=None,
                   help='the jMorp allele-frequency CSV; with it, the cost of NOT '
                        'translating to P groups is measured rather than asserted')
    p.add_argument('--control-group', default='AGP3K',
                   help="the manifest group whose frequencies are compared")
    p.add_argument('--out', required=True)
    return p.parse_args()


def pct(x, dp=2):
    """A share rendered the one way prose is allowed to write it."""
    return f'{100.0 * float(x):.{dp}f} %'


def num(x):
    return f'{int(x):,}'


def main():
    a = parse_args()
    R = Path(a.results)
    f, raw = {}, {}

    def put(key, rendered, value=None):
        f[key] = rendered
        raw[key] = value if value is not None else rendered

    # ---- who was typed -----------------------------------------------------
    man = pd.read_csv(R / '00.manifest/sample_manifest.tsv', sep='\t', dtype=str)
    put('n_typed', num(len(man)), len(man))
    by_group = man['group'].value_counts()
    for g, n in by_group.items():
        put(f'n_typed_{g}', num(n), int(n))
    put('n_platforms', num(man['platform'].nunique()), int(man['platform'].nunique()))

    # ---- the sample chain, all three stages, from the files the .nf reads -----
    # It is three filters and they are NOT the same kind: a workbook flag, a
    # quality decision, an ANCESTRY decision. Collapsing any two of them is the
    # error this block exists to make impossible.
    put('cohort_name', str(a.cohort_name))
    if a.workbook:
        wb = pd.read_excel(a.workbook, engine='openpyxl', dtype=str)
        flag = wb[a.flag_col].astype(str).str.strip().str.lower()
        v6 = wb[flag.isin(('true', '1', 'yes'))]
        put('n_workbook_rows', num(len(wb)), len(wb))
        put('n_v6_flagged', num(len(v6)), len(v6))
        for g, n in v6[a.group_col].value_counts().items():
            put(f'n_v6_{g}', num(n), int(n))

        if a.keep_id:
            qc = {ln.split()[0] for ln in
                  Path(a.keep_id).read_text().splitlines() if ln.strip()}
            v6_qc = v6[v6[a.id_col].astype(str).isin(qc)]
            put('n_qc_pass', num(len(v6_qc)), len(v6_qc))
            put('n_dropped_sample_qc', num(len(v6) - len(v6_qc)), len(v6) - len(v6_qc))
            for g, n in v6_qc[a.group_col].value_counts().items():
                put(f'n_qc_pass_{g}', num(n), int(n))
            # by platform, because the exclusions are inherited and a reader will
            # ask which platforms lost samples.
            for p_, n in v6[~v6[a.id_col].astype(str).isin(qc)][
                    a.platform_col].value_counts().items():
                put(f'n_dropped_sample_qc_{p_}', num(n), int(n))

            if a.cohort_keep:
                ch = {ln.split()[0] for ln in
                      Path(a.cohort_keep).read_text().splitlines() if ln.strip()}
                put('n_dropped_ancestry', num(len(v6_qc) - len(man)), len(v6_qc) - len(man))
                put('n_cohort_outside_qc', num(len(ch - qc)), len(ch - qc))
                out = v6_qc[~v6_qc[a.id_col].astype(str).isin(ch)]
                for g, n in out[a.group_col].value_counts().items():
                    put(f'n_dropped_ancestry_{g}', num(n), int(n))
                for p_, n in out[a.platform_col].value_counts().items():
                    put(f'n_dropped_ancestry_{p_}', num(n), int(n))

    # ---- containment against each downstream cohort ------------------------
    for spec in a.cohort_pheno:
        name, _, path = spec.partition('=')
        ph = pd.read_csv(path, sep='\t', dtype=str)
        ids = set(ph.iloc[:, 1])
        put(f'n_{name}', num(len(ph)), len(ph))
        put(f'n_{name}_not_typed', num(len(ids - set(man['sample_id']))),
            len(ids - set(man['sample_id'])))
        put(f'n_typed_beyond_{name}', num(len(set(man['sample_id']) - ids)),
            len(set(man['sample_id']) - ids))

    # ---- completeness ------------------------------------------------------
    st = pd.read_csv(R / '03.alleles/typing_status.tsv', sep='\t', dtype=SAMPLE_ID)
    for col in ('called', 'hemizygous', 'not_typed', 'failed', 'ambiguous'):
        put(f'n_{col}_gene_calls', num(st[col].sum()), int(st[col].sum()))
    put('n_genes_typed', num(st['genes'].iloc[0]), int(st['genes'].iloc[0]))
    put('n_gene_calls', num(len(st) * st['genes'].iloc[0]), int(len(st) * st['genes'].iloc[0]))

    gq = pd.read_csv(R / '05.qc/typing_qc_gene.tsv', sep='\t')
    put('n_genes_scored', num(len(gq)), len(gq))
    put('n_genes_failed', num((gq['failed'] > 0).sum()), int((gq['failed'] > 0).sum()))
    for _, r in gq.iterrows():
        put(f"call_rate_{r['gene']}", f"{float(r['call_rate']):.4f}", float(r['call_rate']))

    # ---- what actually ran, read back from the run record --------------------
    rm = R / '_run_info/run_manifest.json'
    if rm.is_file():
        man_json = json.loads(rm.read_text())
        if 'hlahd_version' in man_json:
            put('hlahd_version', str(man_json['hlahd_version']))
        # The IMGT release is IN the pinned dictionary path and nowhere else, so
        # it is parsed from there rather than typed beside it.
        m_ = re.search(r'IMGT-HLA-([0-9.]+?)_', str(man_json.get('hlahd_dictionary', '')))
        if m_:
            put('imgt_release', m_.group(1))

    sq = pd.read_csv(R / '05.qc/typing_qc_sample.tsv', sep='\t', dtype=SAMPLE_ID)
    put('n_call_rate_values', num(sq['call_rate'].nunique()), int(sq['call_rate'].nunique()))
    top2 = float(sq['call_rate'].value_counts(normalize=True).head(2).sum())
    put('share_top_two_call_rates', pct(top2, 0), top2)
    put('max_field_depth_seen', num(4 if sq['field4'].sum() else 3),
        int(4 if sq['field4'].sum() else 3))
    put('n_field4_calls', num(sq['field4'].sum()), int(sq['field4'].sum()))

    pq = pd.read_csv(R / '05.qc/typing_qc_platform.tsv', sep='\t')
    # The 2-field share per platform, which the report quotes as a range. Derived
    # the same way plot_typing_qc.py derives it, from the per-sample counts.
    fcols = [c for c in ('field2', 'field3', 'field4') if c in sq.columns]
    if fcols:
        cmp_ = sq.groupby('platform')[fcols].sum()
        f2 = cmp_['field2'] / cmp_[fcols].sum(axis=1)
        put('field2_share_min', pct(f2.min(), 1), float(f2.min()))
        put('field2_share_max', pct(f2.max(), 1), float(f2.max()))
        put('platform_best_resolution', str(f2.idxmin()))
        put('platform_worst_resolution', str(f2.idxmax()))
    for _, r in pq.iterrows():
        k = str(r['platform']).replace(' ', '_')
        put(f'n_{k}', num(r['n_samples']), int(r['n_samples']))
        put(f'call_rate_{k}', f"{float(r['call_rate_mean']):.4f}", float(r['call_rate_mean']))
        if 'mean_field_depth_mean' in pq.columns:
            put(f'field_depth_{k}', f"{float(r['mean_field_depth_mean']):.4f}",
                float(r['mean_field_depth_mean']))

    ad = pd.read_csv(R / '04.residues/allele_dosage.tsv', sep='\t', nrows=1)
    rd = pd.read_csv(R / '04.residues/residue_dosage.tsv', sep='\t', nrows=1)
    dp = pd.read_csv(R / '04.residues/residue_diplotype.tsv', sep='\t', nrows=1)
    put('n_allele_columns', num(len(ad.columns) - 1), len(ad.columns) - 1)
    put('n_residue_columns', num(len(rd.columns) - 1), len(rd.columns) - 1)
    put('n_diplotype_columns', num(len(dp.columns) - 1), len(dp.columns) - 1)

    sites = pd.read_csv(R / '04.residues/residue_sites.tsv', sep='\t')
    put('n_genes_with_residues', num(len(sites)), len(sites))
    put('n_reference_positions', num(sites['positions'].sum()), int(sites['positions'].sum()))
    put('n_polymorphic_positions', num(sites['polymorphic'].sum()),
        int(sites['polymorphic'].sum()))

    # ---- the frequency check, per reference panel --------------------------
    for tag, stem in (('jmorp', 'allele_frequency_summary.tsv'),
                      ('1kg', 'allele_frequency_summary.1kg.tsv')):
        path = R / '05.qc' / stem
        if not path.is_file():
            continue
        s = pd.read_csv(path, sep='\t')
        put(f'n_loci_{tag}', num(len(s)), len(s))
        # A locus whose two denominators mean different things is reported but not
        # counted -- allele_freq_check.py decides that, this only reads it.
        cmpb = (s['comparable'] == 1) if 'comparable' in s.columns else s['gene'].notna()
        put(f'n_loci_comparable_{tag}', num(cmpb.sum()), int(cmpb.sum()))
        put(f'n_loci_not_comparable_{tag}', num((~cmpb).sum()), int((~cmpb).sum()))
        if (~cmpb).any():
            put(f'locus_not_comparable_{tag}', ', '.join(sorted(s.loc[~cmpb, 'gene'])))
            nc = s[~cmpb].iloc[0]
            put(f'frac_chr_ctrl_not_comparable_{tag}', pct(nc['frac_chr_ctrl'], 1),
                float(nc['frac_chr_ctrl']))
            put(f'frac_chr_ref_not_comparable_{tag}', pct(nc['frac_chr_ref'], 0),
                float(nc['frac_chr_ref']))
            # The r the excluded locus would have reported. It is quoted only to
            # be disowned -- it measures the denominator mismatch, not the typing.
            put(f'r_not_comparable_{tag}', f"{float(nc['pearson_r']):.3f}",
                float(nc['pearson_r']))
        # Both denominators at every locus, comparable or not: this pair is what
        # the comparability test reads, so the prose that explains the test has
        # to be able to quote it rather than restate it.
        for _, r in s.iterrows():
            put(f"frac_chr_ctrl_{r['gene']}_{tag}", pct(r['frac_chr_ctrl'], 1),
                float(r['frac_chr_ctrl']))
            put(f"frac_chr_ref_{r['gene']}_{tag}", pct(r['frac_chr_ref'], 1),
                float(r['frac_chr_ref']))
        s = s[cmpb]
        put(f'r_min_{tag}', f"{s['pearson_r'].min():.4f}", float(s['pearson_r'].min()))
        put(f'r_max_{tag}', f"{s['pearson_r'].max():.4f}", float(s['pearson_r'].max()))
        ok = s[s['pearson_r'] >= R_CONCORDANT]
        put('r_concordant_threshold', f'{R_CONCORDANT:.2f}', R_CONCORDANT)
        put(f'n_loci_concordant_{tag}', num(len(ok)), len(ok))
        put(f'r_min_concordant_{tag}', f"{ok['pearson_r'].min():.3f}",
            float(ok['pearson_r'].min()))
        worst = s.loc[s['pearson_r'].idxmin()]
        put(f'worst_locus_{tag}', str(worst['gene']))
        put(f'worst_locus_r_{tag}', f"{float(worst['pearson_r']):.2f}",
            float(worst['pearson_r']))
        put(f'n_chr_ref_{tag}', num(s['n_chr_ref'].max()), int(s['n_chr_ref'].max()))
        put(f'n_ref_individuals_{tag}', num(s['n_chr_ref'].max() // 2),
            int(s['n_chr_ref'].max() // 2))
        for _, r in s.iterrows():
            put(f"r_{r['gene']}_{tag}", f"{float(r['pearson_r']):.4f}", float(r['pearson_r']))
            if 'n_shared' in s.columns:
                put(f"n_shared_{r['gene']}_{tag}", num(r['n_shared']), int(r['n_shared']))
        if 'n_shared' in s.columns:
            put(f'n_shared_min_{tag}', num(s['n_shared'].min()), int(s['n_shared'].min()))
            put(f'n_shared_max_{tag}', num(s['n_shared'].max()), int(s['n_shared'].max()))
            put(f'locus_fewest_shared_{tag}',
                str(s.loc[s['n_shared'].idxmin(), 'gene']))

    # ---- the confound ------------------------------------------------------
    cf = R / '05.qc/typing_confound.tsv'
    if cf.is_file():
        c = pd.read_csv(cf, sep='\t')
        for _, r in c.iterrows():
            k = str(r['stratum']).replace(' ', '_').replace(':', '_').replace('-', '_')
            put(f'unconfirmed_{k}', pct(r['share_unconfirmed']), float(r['share_unconfirmed']))
            put(f'median_depth_{k}', f"{float(r['median_depth']):.2f}", float(r['median_depth']))
        fp = c['fisher_p_matched'].dropna()
        if len(fp):
            put('fisher_p_matched', f'{float(fp.iloc[-1]):.3g}', float(fp.iloc[-1]))
        mw = c['mannwhitney_p_depth'].dropna()
        if len(mw):
            put('mannwhitney_p_depth', f'{float(mw.iloc[0]):.3g}', float(mw.iloc[0]))
        # The group contrast's own P, per panel. The primary panel's is on the
        # `cases` row; each secondary panel adds its own row pair at the end.
        if 'fisher_p_group' in c.columns:
            g = c[['stratum', 'fisher_p_group']].dropna()
            g = g[g.fisher_p_group.astype(str).str.strip() != '']
            for i, (_, r) in enumerate(g.iterrows()):
                key = 'fisher_p_group' if r.stratum == 'cases' else \
                      'fisher_p_group_' + str(r.stratum).split(':')[0] \
                      .replace(' ', '_').replace('-', '_')
                put(key, f'{float(r.fisher_p_group):.3g}', float(r.fisher_p_group))

    fs = R / '05.qc/allele_frequency_sample.tsv'
    if fs.is_file():
        smp = pd.read_csv(fs, sep='\t', dtype=SAMPLE_ID)
        put('n_chr_per_sample_max', num(smp['n_chr'].max()), int(smp['n_chr'].max()))
        put('n_chr_per_sample_min', num(smp['n_chr'].min()), int(smp['n_chr'].min()))
        put('unconfirmed_cohort', pct(smp['n_unconfirmed'].sum() / smp['n_chr'].sum()),
            float(smp['n_unconfirmed'].sum() / smp['n_chr'].sum()))

    # The same contrast under each panel. The absolute level is a property of the
    # panel -- it doubles between them -- while the ratio between the groups
    # barely moves, which is the whole argument for reading only the contrast.
    for tag, stem in (('jmorp', 'allele_frequency_sample.tsv'),
                      ('1kg', 'allele_frequency_sample.1kg.tsv')):
        path = R / '05.qc' / stem
        if not path.is_file():
            continue
        d = pd.read_csv(path, sep='\t', dtype=SAMPLE_ID)
        g = d.groupby('group').agg(chr=('n_chr', 'sum'), unc=('n_unconfirmed', 'sum'))
        g['share'] = g.unc / g.chr
        for grp, r in g.iterrows():
            put(f'unconfirmed_{grp}_{tag}', pct(r.share), float(r.share))
        put(f'unconfirmed_cohort_{tag}', pct(g.unc.sum() / g.chr.sum()),
            float(g.unc.sum() / g.chr.sum()))
        if len(g) == 2:
            hi, lo = g['share'].max(), g['share'].min()
            put(f'unconfirmed_ratio_{tag}', f'{hi / lo:.2f}', float(hi / lo))
        put(f'n_chr_scored_{tag}', num(g.chr.sum()), int(g.chr.sum()))
        for grp, r in g.iterrows():
            put(f'n_chr_{grp}_{tag}', num(r.chr), int(r.chr))
            put(f'n_unconfirmed_{grp}_{tag}', num(r.unc), int(r.unc))
        put(f'n_chr_per_sample_max_{tag}', num(d.n_chr.max()), int(d.n_chr.max()))
        put(f'n_chr_per_sample_min_{tag}', num(d.n_chr.min()), int(d.n_chr.min()))

    # ---- what a non-comparable locus looks like ---------------------------
    # OPEN_QUESTIONS section 3 quotes these, so they are derived rather than
    # measured once and pasted. For each locus the denominator test rejected:
    # the allele that absorbs the mismatch, and the correlation once it is
    # removed and both sides renormalised.
    smy = R / '05.qc/allele_frequency_summary.tsv'
    ckp = R / '05.qc/allele_frequency_check.tsv'
    if smy.is_file() and ckp.is_file():
        sm_ = pd.read_csv(smy, sep='\t')
        ck_ = pd.read_csv(ckp, sep='\t')
        if 'comparable' in sm_.columns:
            for g in sm_.loc[sm_.comparable == 0, 'gene']:
                d_ = ck_[ck_.gene == g].copy()
                d_['adiff'] = (d_.freq_ctrl.astype(float)
                               - d_.freq_ref.astype(float)).abs()
                top = d_.loc[d_.adiff.idxmax()]
                put(f'absorbing_allele_{g}', str(top.allele))
                put(f'absorbing_freq_ref_{g}', f'{float(top.freq_ref):.3f}',
                    float(top.freq_ref))
                put(f'absorbing_freq_ctrl_{g}', f'{float(top.freq_ctrl):.3f}',
                    float(top.freq_ctrl))
                rest = d_[d_.allele != top.allele]
                kc_, kr_ = rest.count_ctrl.sum(), rest.count_ref.sum()
                if kc_ and kr_ and len(rest) > 2:
                    rr = float(np.corrcoef(rest.count_ctrl / kc_,
                                           rest.count_ref / kr_)[0, 1])
                    put(f'r_without_absorbing_{g}', f'{rr:.3f}', rr)

    # ---- the per-locus view of the same contrast --------------------------
    # OPEN_QUESTIONS section 3 is a per-locus table, and it was hand-typed. It is
    # derived here so it cannot outlive the run that produced it.
    ck = R / '05.qc/allele_frequency_check.tsv'
    if ck.is_file():
        from scipy import stats as _st
        chk = pd.read_csv(ck, sep='\t')
        for g, dd in chk.groupby('gene'):
            nc = float(dd.n_chr_ctrl.max() or 0)
            nk = float(dd.n_chr_case.max() or 0)
            uc = float(dd.loc[dd.in_obs_only == 1, 'count_ctrl'].sum())
            uk = float(dd.loc[dd.in_obs_only == 1, 'count_case'].sum())
            if nc:
                put(f'unconfirmed_ctrl_{g}', pct(uc / nc), uc / nc)
                put(f'n_tail_alleles_{g}', num((dd.in_obs_only == 1).sum()),
                    int((dd.in_obs_only == 1).sum()))
                put(f'n_tail_chr_{g}', num(uc), int(uc))
                common = dd[(dd.freq_ref.astype(float) >= 0.05) & (dd.in_obs_only == 0)]
                sf = float((common.freq_ref.astype(float)
                            - common.freq_ctrl.astype(float)).sum())
                put(f'shortfall_{g}', f'{sf:.4f}', sf)
            if nc and nk:
                put(f'confirmed_ctrl_{g}', pct(1 - uc / nc), 1 - uc / nc)
                put(f'confirmed_case_{g}', pct(1 - uk / nk), 1 - uk / nk)
                pv = float(_st.fisher_exact([[uc, nc - uc], [uk, nk - uk]])[1])
                put(f'fisher_p_{g}', f'{pv:.3g}', pv)
                # The DIFFERENCE is what the per-locus table in OPEN_QUESTIONS
                # section 3 actually prints, so it is derived rather than
                # subtracted by hand in the prose.
                d_ = (1 - uk / nk) - (1 - uc / nc)
                put(f'confirmed_diff_{g}', f'{"+" if d_ >= 0 else "-"}{abs(d_) * 100:.2f} %',
                    float(d_))
        # The POOLED version comes from the per-sample table, not from summing the
        # per-gene denominators: n_chr_case is not the same at every locus, and
        # summing per-gene maxima gave 97.88 % where the chromosome-weighted
        # figure the rest of this component quotes is 96.47 %. One denominator,
        # and it is the one typing_confound.tsv uses.
        fs_ = R / '05.qc/allele_frequency_sample.tsv'
        if fs_.is_file():
            d_ = pd.read_csv(fs_, sep='\t', dtype=SAMPLE_ID)
            gg = d_.groupby('group').agg(chr=('n_chr', 'sum'), unc=('n_unconfirmed', 'sum'))
            for grp, r_ in gg.iterrows():
                put(f'confirmed_{grp}_pooled', pct(1 - r_.unc / r_.chr),
                    float(1 - r_.unc / r_.chr))
            if len(gg) == 2:
                lo_, hi_ = sorted(1 - gg.unc / gg.chr)
                put('confirmed_diff_pooled', f'+{(hi_ - lo_) * 100:.2f} %', float(hi_ - lo_))

    # ---- things prose quotes that no other block derives -------------------
    man_ids = man['sample_id'].astype(str)
    put('n_leading_zero_ids', num((man_ids.str.startswith('0')).sum()),
        int((man_ids.str.startswith('0')).sum()))
    put('n_after_numeric_join', num(len(man) - (man_ids.str.startswith('0')).sum()),
        int(len(man) - (man_ids.str.startswith('0')).sum()))
    cf2 = R / '05.qc/typing_confound.tsv'
    if cf2.is_file():
        cc = pd.read_csv(cf2, sep='\t')
        w = cc[cc['stratum'].astype(str).str.startswith('matched')]
        if len(w):
            put('n_matched_window', num(w['n'].sum()), int(w['n'].sum()))

    # The superseded 3,569-sample contrast, so a document can say what the cohort
    # restriction did to it without hand-typing the old number.
    A = R / '_superseded.3569/05.qc'
    for tag, stem in (('jmorp', 'allele_frequency_sample.tsv'),
                      ('1kg', 'allele_frequency_sample.1kg.tsv')):
        path = A / stem
        if not path.is_file():
            continue
        d0 = pd.read_csv(path, sep='\t', dtype=SAMPLE_ID)
        g0 = d0.groupby('group').agg(chr=('n_chr', 'sum'), unc=('n_unconfirmed', 'sum'))
        g0['share'] = g0.unc / g0.chr
        for grp, r in g0.iterrows():
            put(f'superseded_unconfirmed_{grp}_{tag}', pct(r.share), float(r.share))
        if len(g0) == 2:
            hi, lo = g0['share'].max(), g0['share'].min()
            put(f'superseded_unconfirmed_ratio_{tag}', f'{hi / lo:.2f}', float(hi / lo))
        put(f'n_superseded_{tag}', num(len(d0)), len(d0))

    # ---- the residue matrix's own shape, which OUTPUTS.md describes ---------
    # These moved when the cohort did -- a position that was polymorphic in 3,569
    # samples need not be in 3,101 -- so they are derived rather than remembered.
    rdo = pd.read_csv(R / '04.residues/residue_dosage.tsv', sep='\t',
                      dtype=SAMPLE_ID).set_index('sample_id')
    rnum = rdo.apply(pd.to_numeric, errors='coerce')
    DROP = {'DRB3', 'DRB4', 'DRB5'}          # OPEN_QUESTIONS section 1
    by_pos = {}
    for col in rnum.columns:
        gene, position, _aa = col.split(':', 2)
        if gene in DROP:
            continue
        by_pos.setdefault(f'{gene}:{position}', []).append(col)
    psum = pd.DataFrame({k: rnum[v].sum(axis=1) for k, v in by_pos.items()})
    eq2 = (psum == 2)
    und = 1 - eq2.mean(axis=0)
    put('n_usable_positions', num(len(by_pos)), len(by_pos))
    put('share_positions_summing_to_2', pct(eq2.to_numpy().mean()),
        float(eq2.to_numpy().mean()))
    put('median_undetermined_rate', pct(und.median(), 4), float(und.median()))
    put('n_positions_over_5pct', num((und > 0.05).sum()), int((und > 0.05).sum()))
    put('worst_undetermined_position', str(und.idxmax()))
    put('worst_undetermined_rate', pct(und.max()), float(und.max()))

    pg = R / '05.qc/allele_pgroup_map.tsv'
    if pg.is_file():
        m = pd.read_csv(pg, sep='\t')
        put('n_alleles_carried', num(len(m)), len(m))
        put('n_alleles_not_in_reference', num((m['in_reference'] == 0).sum()),
            int((m['in_reference'] == 0).sum()))
        # What collapsing allele_dosage.tsv onto P groups would cost, which
        # METHODS section 13 quotes as the reason not to.
        put('n_pgroups_carried', num(m['p_group'].nunique()), int(m['p_group'].nunique()))
        # What the harmonisation actually costs, at the loci that are compared.
        # An allele whose name is unchanged met the reference at the full 2-field
        # name; one that changed was merged, and only where the reference had
        # merged it too.
        LOCI = set(pd.read_csv(R / '05.qc/allele_frequency_summary.tsv',
                               sep='\t')['gene'])
        mm = m[m.gene.isin(LOCI)]
        unchanged = int((mm.allele == mm.p_group).sum())
        put('n_alleles_at_compared_loci', num(len(mm)), len(mm))
        put('n_alleles_name_unchanged', num(unchanged), unchanged)
        put('n_alleles_merged', num(len(mm) - unchanged), len(mm) - unchanged)
        grp = mm.groupby('p_group').allele.nunique()
        coll = grp[grp > 1]
        put('n_alleles_in_multi_groups', num(coll.sum()), int(coll.sum()))
        put('n_multi_groups', num(len(coll)), len(coll))
        put('n_names_lost_to_merging', num(int(coll.sum()) - len(coll)),
            int(coll.sum()) - len(coll))

        # The counterfactual behind "the translation is not optional": join our
        # raw two-field names straight against the reference's own names and see
        # what is left. The reference publishes a MIXTURE of P-group and plain
        # names, so a raw call joins only where the reference happens to write
        # that allele plain -- which is why this is measured and not guessed.
        if a.reference_jmorp:
            ref_names = set(pd.read_csv(a.reference_jmorp)['allele'].astype(str))
            # The mixed convention itself, counted rather than described: the
            # report opens the harmonisation section with this split, so it has
            # to move if the panel is ever replaced.
            n_p = sum(1 for n in ref_names if n.endswith('P'))
            put('n_ref_names_jmorp', num(len(ref_names)), len(ref_names))
            put('n_ref_names_pgroup_jmorp', num(n_p), n_p)
            put('n_ref_names_plain_jmorp', num(len(ref_names) - n_p),
                len(ref_names) - n_p)
            # How many of OUR alleles the report's worked example collapses,
            # so the sentence explaining a P group counts rather than asserts.
            worked = f.get('absorbing_allele_DRB3')
            if worked:
                n_in = int((m.p_group == worked).sum())
                put('n_alleles_in_worked_group', num(n_in), n_in)
            ctrl = set(man.loc[man['group'] == a.control_group, 'sample_id'])
            dos = pd.read_csv(R / '04.residues/allele_dosage.tsv', sep='\t',
                              dtype=SAMPLE_ID)
            chrom = dos[dos.sample_id.isin(ctrl)].drop(columns=['sample_id']).sum()
            gene_of = dict(zip(m.allele, m.gene))
            cols = [c for c in chrom.index if gene_of.get(c) in LOCI]
            tot = float(chrom[cols].sum())
            hit = float(chrom[[c for c in cols if c in ref_names]].sum())
            put('unmatched_without_pgroups_jmorp', pct(1.0 - hit / tot),
                1.0 - hit / tot)
            dead = []
            for g in sorted(LOCI):
                gc = [c for c in cols if gene_of[c] == g]
                if float(chrom[gc].sum()) and not float(
                        chrom[[c for c in gc if c in ref_names]].sum()):
                    dead.append(g)
            put('n_loci_unmatched_without_pgroups_jmorp', num(len(dead)), len(dead))
            put('loci_unmatched_without_pgroups_jmorp', ', '.join(dead))

    out = {'_note': 'Generated by build_facts.py. Every number this component writes into prose '
                    'comes from here; nothing is hand-typed. See verify.sh.',
           'rendered': f, 'raw': raw}
    Path(a.out).write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(f'[build_facts] {len(f)} fact(s) -> {a.out}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in build_facts: {e}', file=sys.stderr)
        sys.exit(1)
