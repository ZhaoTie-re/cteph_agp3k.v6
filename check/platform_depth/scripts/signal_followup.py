#!/usr/bin/env python3
"""Follow-up of what the checks leave open: counts and frequencies, no new association test (report stage).

Every table holds the three cohorts (column `cohort`; the cohorts are nested, so each is a selection of the
full-mainland samples). Groups: each case platform (and, for the GWAS leads, its registration study), each
control source; the pooled groups beside them.

  followup.gwas_af.tsv          every GWAS lead: frequency of the variant's ALT allele per group, with ToMMo
                                61KJPN and 1KG JPT (integration/saige_signals)
  followup.gene_carriers.tsv    every tested variant of the Tier 1-2 genes (both sets) and each gene as a whole:
                                carriers per group; the genotypes the cohort's filters set missing although the
                                pre-QC record carries the ALT (non-reference GT or >= 2 ALT reads: possible
                                carriers hidden as missing); ToMMo, gnomAD EAS, 1KG EAS / JPT frequencies
                                (integration/gene_signals) and the carriers they predict
  followup.gene_filters.tsv     the cohort's genotype filters re-applied to the pre-QC records against the published
                                genotypes, per variant (must agree for the audit to stand)
  followup.hla_v57.tsv          every down-sampled case whose DQB1:57:V dosage differs between any two arms: dosage
                                per arm, the DQB1 call, its read support and HLA-HD confidence
  followup.lead_preqc.tsv       every GWAS lead x group: genotypes before and after the cohort's filters (pre-QC
                                carriers, missing, set missing by the filters, set missing although the record
                                carries the ALT): can missingness hide carriers at a lead?
  followup.lead_window.tsv      every pre-QC record within 50 bp of a GWAS lead: carriers per group (another
                                representation of the lead's allele, or a linked variant, a platform may call)
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402

# the cohort's genotype filters (wgs.auto.par/select.auto.par.v6.nf:252-260)
MIN_GQ, MIN_DP, HET_AF = 20, 8, (0.2, 0.8)


def groups(results, ids):
    """per sample: case?, the platform / control-source group, and the case groups split by registration study"""
    st = pd.read_csv(results / '00.samples' / 'sample_table.tsv', sep='\t', dtype=str)
    info = pd.read_excel(K.INFO_XLSX, usecols=['ID_JHRPv6', 'Project', 'Registration_Study'])
    info = info.dropna(subset=['ID_JHRPv6']).drop_duplicates('ID_JHRPv6')
    x = st.merge(info, left_on='sample_id', right_on='ID_JHRPv6', how='left').set_index('sample_id').reindex(ids)
    case = x.status == 'case'
    plat = np.where(case, RL.CASE + x.platform.fillna('NA'), RL.CONTROL + x.Project.map(RL.CONTROL_SOURCES).fillna('NA'))
    study = np.where(case, plat + ' | ' + x.Registration_Study.map(RL.STUDY).fillna('NA') + ' study', plat)
    return case.values, plat, study


def gwas_af(a, cohort):
    ids, vids, dose = RL.alt_dosage(a.results / '04.genotypes' / 'original' / cohort / 'gwas')
    case, plat, study = groups(a.results, ids)
    leads = RL.gwas_leads(RL.read_signals(a.results))
    ref = pd.read_csv(a.saige_signals, sep='\t', dtype={'variant_id': str}, low_memory=False)
    ref = ref.drop_duplicates('variant_id').set_index('variant_id')
    rows = []
    for u in leads.itertuples():
        j = int(np.where(vids == u.variant_id)[0][0])
        d = dose[j]
        sets = [(RL.CASE_ALL, case), (RL.CONTROL_ALL, ~case)]
        sets += [(g, plat == g) for g in sorted(set(plat))] + [(g, study == g) for g in sorted(set(study) - set(plat))]
        for g, m in sets:
            ok = m & np.isfinite(d)
            rows.append({'cohort': cohort, 'signal': u.unit, 'variant_id': u.variant_id, 'group': g,
                         'n': int(ok.sum()), 'alt_freq': d[ok].mean() / 2 if ok.any() else np.nan})
        for col, lab in [('tommo_af', RL.REF_TOMMO), ('kg_jpt_af', RL.REF_JPT)]:
            v = ref.loc[u.variant_id, col] if u.variant_id in ref.index and col in ref.columns else np.nan
            rows.append({'cohort': cohort, 'signal': u.unit, 'variant_id': u.variant_id, 'group': lab, 'n': np.nan,
                         'alt_freq': pd.to_numeric(v, errors='coerce')})
    return pd.DataFrame(rows)


def set_positions(genes, cohort):
    """variant positions of each gene in each set, from the cohort's published rvtest set files"""
    out = {}
    for st in RL.VARIANT_SETS:
        for f in sorted((K.ASSOC_GENE / cohort / '02.test_inputs' / st).glob('rvtest.set.chr*.txt')):
            for ln in f.read_text().splitlines():
                g, regs = ln.split('\t')[:2]
                if g in genes:
                    out[(g, st)] = {(f'chr{r.split(":")[0]}', int(r.split(':')[1].split('-')[0])) for r in regs.split(',')}
    return out


def pre_qc(a, vids):
    """GT, AD, DP, GQ and FORMAT/AF of every sample at the variants, from the pre-QC cohort VCF"""
    recs = []
    for v in vids:
        ch, pos, ref, alt = v.split(':')
        vcf = a.cohort_vcf.replace('@CHROM@', ch)
        q = subprocess.run([a.bcftools, 'query', '-r', f'{ch}:{pos}',
                            '-f', r'[%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%AD\t%DP\t%GQ\t%AF\n]', vcf],
                           capture_output=True, text=True, check=True)
        n0 = len(recs)
        for ln in q.stdout.splitlines():
            p, r_, al, s, gt, ad, dp, gq, af = ln.split('\t')
            if p == pos and r_ == ref and al == alt:
                recs.append((v, s, gt, ad, dp, gq, af))
        if len(recs) == n0:
            sys.exit(f'ABORT: {v}: no biallelic record in {vcf}')
    d = pd.DataFrame(recs, columns=['variant_id', 'sample_id', 'GT', 'AD', 'DP', 'GQ', 'AFf'])
    num = lambda s: pd.to_numeric(s.replace('.', np.nan), errors='coerce')  # noqa: E731
    gt = d['GT'].str.replace('|', '/', regex=False)
    pre = np.select([gt == '0/0', gt.isin(['0/1', '1/0']), gt == '1/1'], [0.0, 1.0, 2.0], np.nan)
    dp, gq = num(d['DP']), num(d['GQ'])
    af = pd.to_numeric(d['AFf'].str.split(',').str[0], errors='coerce')
    fail = (gq < MIN_GQ) | (dp < MIN_DP) | ((pre == 1) & ((af < HET_AF[0]) | (af > HET_AF[1]))) | ((pre > 0) & af.isna())
    d['pre'] = pre
    d['post'] = np.where(fail, np.nan, pre)
    d['alt_reads'] = d['AD'].str.split(',').apply(lambda x: sum(int(v) for v in x[1:] if v.isdigit()))
    d['hidden'] = np.isnan(d['post']) & ((d['pre'] > 0) | (d['alt_reads'] >= 2))
    return d


WINDOW = 50


def lead_preqc(a, cohort, pre_cache):
    """per lead and group: the pre-QC genotypes against the published ones (the gene audit applied to the leads)"""
    ids, vids, dose = RL.alt_dosage(a.results / '04.genotypes' / 'original' / cohort / 'gwas')
    case, plat, _ = groups(a.results, ids)
    leads = RL.gwas_leads(RL.read_signals(a.results))
    new = [v for v in leads.variant_id if v not in pre_cache]
    if new:
        pre_cache.update(dict(tuple(pre_qc(a, new).groupby('variant_id'))))
    sets = [(RL.CASE_ALL, case), (RL.CONTROL_ALL, ~case)] + [(g, plat == g) for g in sorted(set(plat))]
    rows = []
    for u in leads.itertuples():
        q = pre_cache[u.variant_id].set_index('sample_id').reindex(ids)
        pub = dose[int(np.where(vids == u.variant_id)[0][0])]
        pre = q['pre'].to_numpy(float)
        alt = q['alt_reads'].fillna(0).to_numpy(float)
        for g, m in sets:
            rows.append({'cohort': cohort, 'signal': u.unit, 'variant_id': u.variant_id, 'group': g, 'n': int(m.sum()),
                         'pre_carriers': int((m & (pre > 0)).sum()), 'pre_missing': int((m & np.isnan(pre)).sum()),
                         'published_carriers': int((m & (pub > 0)).sum()),
                         'published_missing': int((m & np.isnan(pub)).sum()),
                         'set_missing_by_filters': int((m & ~np.isnan(pre) & np.isnan(pub)).sum()),
                         'set_missing_with_alt': int((m & np.isnan(pub) & ((pre > 0) | (alt >= 2))).sum()),
                         'alt_read_samples': int((m & (alt >= 1)).sum())})
    return pd.DataFrame(rows)


def lead_window(a, cohort):
    """every pre-QC record within WINDOW bp of a lead: carriers (any non-reference GT) per group"""
    ids, _, _ = RL.alt_dosage(a.results / '04.genotypes' / 'original' / cohort / 'gwas')
    case, plat, _ = groups(a.results, ids)
    col = {s: i for i, s in enumerate(ids)}
    sets = [(RL.CASE_ALL, case), (RL.CONTROL_ALL, ~case)] + [(g, plat == g) for g in sorted(set(plat))]
    rows = []
    for u in RL.gwas_leads(RL.read_signals(a.results)).itertuples():
        ch, pos = u.variant_id.split(':')[:2]
        q = subprocess.run([a.bcftools, 'query', '-r', f'{ch}:{int(pos) - WINDOW}-{int(pos) + WINDOW}',
                            '-f', r'%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n', a.cohort_vcf.replace('@CHROM@', ch)],
                           capture_output=True, text=True, check=True)
        for ln in q.stdout.splitlines():
            f = ln.split('\t')
            car = np.zeros(len(ids), bool)
            for x in f[3:]:
                smp, gt = x.split('=')
                if smp in col and any(al not in ('0', '.') for al in gt.replace('|', '/').split('/')):
                    car[col[smp]] = True
            vid = f'{ch}:{f[0]}:{f[1]}:{f[2]}'
            for g, m in sets:
                rows.append({'cohort': cohort, 'signal': u.unit, 'lead': u.variant_id, 'variant_id': vid,
                             'is_lead': vid == u.variant_id, 'offset_bp': int(f[0]) - int(pos), 'group': g,
                             'n': int(m.sum()), 'carriers': int((car & m).sum())})
    return pd.DataFrame(rows)


def gene_carriers(a, cohort, pre_cache):
    ids, vids, dose = RL.alt_dosage(a.results / '04.genotypes' / 'original' / cohort / 'gene')
    case, plat, _ = groups(a.results, ids)
    rv = pd.read_csv(a.results / '00.samples' / 'sample_table.tsv', sep='\t', dtype=str)
    keep = np.isin(ids, rv.loc[rv[f'rv_{cohort}'] == '1', 'sample_id'])
    sig = RL.read_signals(a.results)
    genes = sorted(sig.loc[sig.kind == 'gene', 'unit'].unique())
    sets = set_positions(genes, cohort)
    pos = {(v.split(':')[0], int(v.split(':')[1])): i for i, v in enumerate(vids)}
    want = sorted({vids[pos[p]] for s in sets.values() for p in s if p in pos})
    new = [v for v in want if v not in pre_cache]        # the pre-QC records are the same call set in every cohort
    if new:
        pre_cache.update(dict(tuple(pre_qc(a, new).groupby('variant_id'))))
    pq = pd.concat([pre_cache[v] for v in want]).set_index(['variant_id', 'sample_id'])
    fq = a.gene_freq
    ext = pd.read_csv(fq / 'tommo.tsv', sep='\t', dtype={'variant_id': str}).set_index('variant_id')[['tommo_af']]
    ext = ext.join(pd.read_csv(fq / 'annovar.tsv', sep='\t', dtype={'variant_id': str}).set_index('variant_id')[['gnomad_eas_af']],
                   how='outer')
    ext = ext.join(pd.read_csv(fq / 'kg.tsv', sep='\t', dtype={'variant_id': str}).set_index('variant_id')[['kg_eas_af', 'kg_jpt_af']],
                   how='outer')
    ext = ext.apply(pd.to_numeric, errors='coerce')
    grp = [(RL.CASE_ALL, case), (RL.CONTROL_ALL, ~case)] + [(g, plat == g) for g in sorted(set(plat))]
    rows, filt = [], []
    for (g, st), ps in sorted(sets.items()):
        idx = [pos[p] for p in sorted(ps) if p in pos]
        units = [(vids[i], [i]) for i in idx] + [(f'{g} (gene)', idx)]
        for unit, ii in units:
            car = np.nanmax(np.where(np.isnan(dose[ii]), -1, dose[ii]), axis=0) > 0
            hid = np.zeros(len(ids), bool)
            for i in ii:
                h = pq.loc[vids[i], 'hidden'].reindex(ids).fillna(False).to_numpy(bool)
                hid |= h & ~car
            two_af = sum(2 * ext.loc[vids[i], 'tommo_af'] if vids[i] in ext.index else np.nan for i in ii)
            for gname, m in grp:
                m = m & keep
                n = int(m.sum())
                r = {'cohort': cohort, 'gene': g, 'stratum': st, 'unit': unit, 'group': gname, 'n': n,
                     'carriers': int((car & m).sum()), 'carrier_share': (car & m).sum() / n if n else np.nan,
                     'hidden_as_missing': int((hid & m).sum()), 'tommo_expected_carriers': n * two_af}
                if len(ii) == 1:
                    for c in ext.columns:
                        r[c] = ext.loc[vids[ii[0]], c] if vids[ii[0]] in ext.index else np.nan
                rows.append(r)
    for v in want:
        x = pq.loc[v]['post'].reindex(ids)
        pub = pd.Series(dose[pos[(v.split(':')[0], int(v.split(':')[1]))]], index=ids)
        both = keep & x.notna().to_numpy() & pub.notna().to_numpy()
        filt.append({'cohort': cohort, 'variant_id': v, 'samples': int(keep.sum()),
                     'missing_published': int(pub[keep].isna().sum()),
                     'missing_reapplied': int(x[keep].isna().sum()),
                     'genotype_agreement': float((x[both] == pub[both]).mean()) if both.any() else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(filt)


def hla_v57(a, cohort):
    lab = a.label
    deep = set(pd.read_csv(a.results / '03.calls' / lab / 'deep_cases.tsv', sep='\t', dtype=str).iloc[:, 0])
    col = RL.V57

    def dosage(f):
        x = pd.read_csv(f, sep='\t', dtype=str)
        x = x[x[x.columns[0]].isin(deep)].set_index(x.columns[0])
        return pd.to_numeric(x[col], errors='coerce')
    arms = ['baseline'] + sorted(p.name for p in (a.results / '05.hla' / lab).glob('r[0-9][0-9]'))
    m = pd.DataFrame({ar: dosage(a.results / '05.hla' / lab / ar / 'matrices' / cohort / 'residue_dosage.tsv') for ar in arms})
    m['original'] = dosage(a.hla_v2 / 'results' / '06.assoc_inputs' / cohort / 'residue_dosage.tsv')
    cols = ['original'] + arms
    chg = m[m[cols].apply(lambda r: r.astype(str).nunique() > 1, axis=1)]
    calls = pd.read_csv(a.results / '07.downsampling' / lab / 'hla' / 'hla.calls_long.tsv', sep='\t', dtype=str)
    calls = calls[(calls.gene == 'DQB1') & (calls.comparison == 'replicate_vs_baseline')].set_index(['sample_id', 'arm'])
    base = calls.reset_index().drop_duplicates('sample_id').set_index('sample_id')
    rows = []
    for s, r in chg.iterrows():
        for ar in cols:
            c = calls.loc[(s, ar)] if (s, ar) in calls.index else None
            rows.append({'cohort': cohort, 'sample_id': s, 'arm': ar, 'v57_dosage': r[ar],
                         'change_vs_baseline': r[ar] - r['baseline'],
                         'dqb1_call': (c.alleles_test if c is not None else base.loc[s, 'alleles_ref'] if s in base.index else np.nan),
                         'call_state': c.state_test if c is not None else np.nan,
                         'dqb1_reads': c.reads_test if c is not None else (base.loc[s, 'reads_ref'] if s in base.index else np.nan),
                         'confidence': c.confidence_test if c is not None else np.nan})
    tot = {ar: float(m[ar].sum()) for ar in cols}
    print(f'[signal_followup] {cohort}: DQB1:57:V copies among the down-sampled cases: ' +
          ', '.join(f'{k} {v:.0f}' for k, v in tot.items()))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', type=Path, default=K.ROOT / 'results')
    ap.add_argument('--label', default='downsampling')
    ap.add_argument('--saige-signals', type=Path, default=K.SAIGE_SIGNALS / 'saige_signals.long.tsv')
    ap.add_argument('--gene-freq', type=Path, default=K.PROJECT / 'integration' / 'gene_signals' / 'results' / '03.frequencies')
    ap.add_argument('--hla-v2', type=Path, default=K.PROJECT / 'hla.typing.v2')
    ap.add_argument('--cohort-vcf', default=str(K.PROJECT / 'wgs.auto.par/results/03_annotate_af_norm_gt/@CHROM@.vqc.af.gtnorm.vcf.gz'))
    ap.add_argument('--bcftools', default='bcftools')
    ap.add_argument('--out-dir', type=Path, default=Path('.'))
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    pre_cache = {}
    g = pd.concat([gwas_af(a, c) for c in K.COHORTS], ignore_index=True)
    g.to_csv(a.out_dir / 'followup.gwas_af.tsv', sep='\t', index=False)
    cf = [gene_carriers(a, c, pre_cache) for c in K.COHORTS]
    c = pd.concat([x[0] for x in cf], ignore_index=True)
    f = pd.concat([x[1] for x in cf], ignore_index=True)
    c.to_csv(a.out_dir / 'followup.gene_carriers.tsv', sep='\t', index=False)
    f.to_csv(a.out_dir / 'followup.gene_filters.tsv', sep='\t', index=False)
    bad = f[f.genotype_agreement < 1]
    if len(bad):
        print(f'[signal_followup] WARNING: re-applied filters differ from the published genotypes at {len(bad)} rows')
    h = pd.concat([hla_v57(a, c_) for c_ in K.COHORTS], ignore_index=True)
    h.to_csv(a.out_dir / 'followup.hla_v57.tsv', sep='\t', index=False)
    lp = pd.concat([lead_preqc(a, c_, pre_cache) for c_ in K.COHORTS], ignore_index=True)
    lp.to_csv(a.out_dir / 'followup.lead_preqc.tsv', sep='\t', index=False)
    lw = pd.concat([lead_window(a, c_) for c_ in K.COHORTS], ignore_index=True)
    lw.to_csv(a.out_dir / 'followup.lead_window.tsv', sep='\t', index=False)
    print(f'[signal_followup] three cohorts: {g.signal.nunique()} GWAS leads, {c.unit.nunique()} gene units, '
          f'{h.sample_id.nunique()} down-sampled cases with a DQB1:57:V change; filters agree at '
          f'{int((f.genotype_agreement == 1).sum())}/{len(f)} cohort x variant rows')


if __name__ == '__main__':
    main()
