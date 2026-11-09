#!/usr/bin/env python3
"""Genotype concordance of the re-called arms (user request 2026-10-30: a reported result).

Genotypes are compared at every fixed site (GWAS signals, their regional context, the genes' tested
variants), as ALT-allele counts 0 / 1 / 2 or missing:

  recall_vs_original      baseline (f = 1, new processing) vs the cohort's genotypes: re-processing
  origgvcf_vs_original    the surviving ORIGINAL GVCFs joint-genotyped by our procedure vs the cohort:
                          the joint-genotyping step alone
  recall_vs_origgvcf      baseline vs those original GVCFs: the HaplotypeCaller re-run alone
  replicate_vs_baseline   each replicate vs baseline: the depth effect (same processing, fewer reads)
  replicate_vs_original   each replicate vs the cohort's genotypes: the total change

Per comparison x arm x sample x variant class x MAF class: the 4 x 4 confusion counts. Summaries:
concordance among genotypes called in both, non-reference concordance (1 - NRD, NRD = discordant
non-reference / all with a non-reference call in either), het and hom-alt retention (called the same
in the test), loss to missing, false-het rate, call rates; replicates as median and range.

`reasons`: the share of genotypes the cohort's QC rules no-call, by rule (GQ < 20, DP < 8, het AF
outside 0.2-0.8, AF NaN), for the deep cases in each arm and for the fidelity controls — and, from
the cohort's own pre-QC VCF, for every original case and control group. After down-sampling, the deep
cases should move toward the controls (the positive control that depth is emulated).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bedio as B  # noqa: E402

LAB = {0: '0', 1: '1', 2: '2', -1: 'miss'}
MAF_CLASS = [(0, 0.01, '<1%'), (0.01, 0.05, '1-5%'), (0.05, 0.51, '>=5%')]


def vcf_alt_counts(path, ids):
    """{sample: {variant_id: alt count or -1}} and the FT reasons {sample: {vid: ft}}."""
    vf = pysam.VariantFile(path)
    samples = list(vf.header.samples)
    g = {s: {} for s in samples}
    ft = {s: {} for s in samples}
    for rec in vf.fetch():
        if rec.id not in ids:
            continue
        for s in samples:
            x = rec.samples[s]
            gt = x.get('GT')
            g[s][rec.id] = -1 if gt is None or any(a is None for a in gt) else sum(1 for a in gt if a == 1)
            f = x.get('FT') if 'FT' in rec.format else None
            ft[s][rec.id] = ','.join(f) if isinstance(f, (list, tuple)) else (f or 'PASS')
    return g, ft


def bed_alt_counts(prefix, samples):
    fam = B.read_fam(prefix)
    bim = B.read_bim(prefix)
    codes = B.unpack(B.read_rows(prefix, list(range(len(bim))), len(fam)), len(fam))
    out = {}
    for i, f in enumerate(bim):
        _, _, ref, alt = f[1].split(':', 3)
        c = codes[i]
        if (f[4], f[5]) == (alt, ref):
            alt_n = np.where(c == B.HOM5, 2, np.where(c == B.HET, 1, np.where(c == B.HOM6, 0, -1)))
        elif (f[4], f[5]) == (ref, alt):
            alt_n = np.where(c == B.HOM5, 0, np.where(c == B.HET, 1, np.where(c == B.HOM6, 2, -1)))
        else:
            sys.exit(f'ABORT: {f[1]} has .bim alleles {f[4]}/{f[5]}')
        out[f[1]] = (alt_n, {s: k for k, s in enumerate(fam)})
    return out


def maf_of(orig, vid):
    a, _ = orig[vid]
    called = a[a >= 0]
    if not len(called):
        return np.nan
    p = called.sum() / (2 * len(called))
    return min(p, 1 - p)


def confusion(ref, test):
    m = np.zeros((4, 4), int)
    for r, t in zip(ref, test):
        m[r if r >= 0 else 3, t if t >= 0 else 3] += 1
    return m


def metrics(m):
    """m[ref, test] over codes 0, 1, 2, missing."""
    both = m[:3, :3]
    nb = both.sum()
    disc_nonref = both.sum() - both[0, 0] - both[1, 1] - both[2, 2]
    nonref = both.sum() - both[0, 0]
    return {
        'n': int(m.sum()), 'both_called': int(nb),
        'concordance': both.trace() / nb if nb else np.nan,
        'nonref_concordance': 1 - disc_nonref / nonref if nonref else np.nan,
        'het_retained': m[1, 1] / m[1, :].sum() if m[1, :].sum() else np.nan,
        'homalt_retained': m[2, 2] / m[2, :].sum() if m[2, :].sum() else np.nan,
        'het_lost_to_missing': m[1, 3] / m[1, :].sum() if m[1, :].sum() else np.nan,
        'false_het_rate': m[0, 1] / m[0, :3].sum() if m[0, :3].sum() else np.nan,
        'callrate_ref': m[:3, :].sum() / m.sum() if m.sum() else np.nan,
        'callrate_test': m[:, :3].sum() / m.sum() if m.sum() else np.nan,
    }


def cmd_variants(a):
    fs = pd.read_csv(a.fixed, sep='\t', dtype=str)
    fs['kind'] = fs.kind.str.split(';').map(lambda k: 'gene' if 'gene' in k else k[0])
    ids = set(fs.variant_id)
    plan = pd.read_csv(a.plan, sep='\t', dtype=str)
    calls, _ = vcf_alt_counts(a.vcf, ids)
    orig = {}
    for p in a.original:
        for vid, v in bed_alt_counts(p, None).items():
            orig.setdefault(vid, v)
    missing = ids - set(orig)
    if missing:   # regional-context variants tested only in the narrower cohorts' random-model sets
        print(f'[concordance] {len(missing)} fixed sites are not in the full-mainland filesets; left out')
        if any(k != 'gwas_context' for k in fs[fs.variant_id.isin(missing)].kind):
            sys.exit('ABORT: a signal or gene variant has no original genotypes')
        fs = fs[~fs.variant_id.isin(missing)]
        ids = set(fs.variant_id)
    maf = {v: maf_of(orig, v) for v in ids}
    fs['maf_class'] = [next(lbl for lo, hi, lbl in MAF_CLASS if lo <= maf[v] < hi) for v in fs.variant_id]
    klass = dict(zip(fs.variant_id, zip(fs.kind, fs.maf_class)))

    def orig_vec(sid, vids):
        return np.array([orig[v][0][orig[v][1][sid]] if sid in orig[v][1] else -1 for v in vids])

    def call_vec(name, vids):
        return np.array([calls[name][v] for v in vids])

    plat = dict(zip(plan.sample_id, plan.platform))
    role = dict(zip(plan.sample_id, plan.role))
    arms = sorted(set(plan.arm[plan.role == 'deep']) - {'baseline'})
    vids = sorted(ids)
    rows = []

    def add(comp, arm, sid, ref, test):
        for (kind, mc) in sorted(set(klass[v] for v in vids)):
            sel = np.array([klass[v] == (kind, mc) for v in vids])
            m = confusion(ref[sel], test[sel])
            rows.append({'comparison': comp, 'arm': arm, 'sample_id': sid, 'role': role.get(sid, ''),
                         'platform': plat.get(sid, ''), 'kind': kind, 'maf_class': mc,
                         **{f'n_{LAB[i if i < 3 else -1]}_{LAB[j if j < 3 else -1]}': int(m[i, j])
                            for i in range(4) for j in range(4)}})

    for sid in sorted(set(plan.sample_id)):
        o = orig_vec(sid, vids)
        base = call_vec(f'{sid}__baseline', vids)
        add('recall_vs_original', 'baseline', sid, o, base)
        if f'{sid}__origgvcf' in calls:
            og = call_vec(f'{sid}__origgvcf', vids)
            add('origgvcf_vs_original', 'origgvcf', sid, o, og)
            add('recall_vs_origgvcf', 'baseline', sid, og, base)
        if role.get(sid) == 'deep':
            for arm in arms:
                t = call_vec(f'{sid}__{arm}', vids)
                add('replicate_vs_baseline', arm, sid, base, t)
                add('replicate_vs_original', arm, sid, o, t)
    long = pd.DataFrame(rows)
    long.to_csv(a.out_long, sep='\t', index=False)

    def summarize(d, keys):
        out = []
        cols = [c for c in d.columns if c.startswith('n_')]
        for k, g in d.groupby(keys):
            m = np.zeros((4, 4), int)
            for i, ri in enumerate(['0', '1', '2', 'miss']):
                for j, cj in enumerate(['0', '1', '2', 'miss']):
                    m[i, j] = g[f'n_{ri}_{cj}'].sum()
            out.append({**dict(zip(keys, k if isinstance(k, tuple) else (k,))), **metrics(m),
                        'samples': g.sample_id.nunique()})
        return pd.DataFrame(out)

    by_arm = summarize(long, ['comparison', 'role', 'arm', 'kind', 'maf_class'])
    by_arm_all = summarize(long.assign(kind='all', maf_class='all'), ['comparison', 'role', 'arm', 'kind', 'maf_class'])
    by_plat = summarize(long.assign(kind='all', maf_class='all'), ['comparison', 'role', 'platform', 'arm', 'kind', 'maf_class'])
    by_sample = summarize(long.assign(kind='all', maf_class='all'), ['comparison', 'role', 'platform', 'arm', 'sample_id'])
    tab = pd.concat([by_arm_all, by_arm, by_plat], ignore_index=True)
    tab.to_csv(a.out_summary, sep='\t', index=False)
    by_sample.to_csv(a.out_samples, sep='\t', index=False)
    # replicates: median and range of each metric over the replicate arms
    rep = tab[tab.comparison.str.startswith('replicate_')]
    keys = ['comparison', 'role', 'platform', 'kind', 'maf_class']
    rep = rep.assign(platform=rep.platform.fillna('all'))
    agg = []
    for k, g in rep.groupby(keys, dropna=False):
        row = dict(zip(keys, k))
        for c in ('concordance', 'nonref_concordance', 'het_retained', 'homalt_retained', 'het_lost_to_missing',
                  'false_het_rate', 'callrate_ref', 'callrate_test'):
            v = g[c].dropna()
            row[f'{c}_median'], row[f'{c}_min'], row[f'{c}_max'] = (
                (float(v.median()), float(v.min()), float(v.max())) if len(v) else (np.nan,) * 3)
        row['replicates'] = g.arm.nunique()
        agg.append(row)
    pd.DataFrame(agg).to_csv(a.out_replicates, sep='\t', index=False)
    print(by_arm_all.round(5).to_string(index=False))


def rule_reasons(gt, gq, dp, af):
    """The cohort's genotype rules on CALLED genotypes (a genotype can fail several)."""
    called = gt >= 0
    het = gt == 1
    with np.errstate(invalid='ignore'):
        return {'LowGQ': called & (gq < 20), 'LowDP': called & (dp < 8),
                'ABB_outlier': het & ((af < 0.2) | (af > 0.8)), 'ABB_NaN': called & np.isnan(af)}


def cmd_reasons(a):
    fs = pd.read_csv(a.fixed, sep='\t', dtype=str)
    ids = set(fs.variant_id)
    plan = pd.read_csv(a.plan, sep='\t', dtype=str)
    st = pd.read_csv(a.samples, sep='\t', dtype=str)
    rows = []
    # arms: FT written by VariantFiltration
    _, ft = vcf_alt_counts(a.vcf, ids)
    for name, d in ft.items():
        sid, arm = name.split('__')
        v = pd.Series(d)
        n = len(v)
        r = {'source': 'recalled', 'group': f'{plan.loc[plan.sample_id == sid, "role"].iloc[0]}', 'arm': arm,
             'sample_id': sid, 'n': n}
        for rule in ('LowGQ', 'LowDP', 'ABB_outlier', 'ABB_NaN'):
            r[rule] = float(v.str.contains(rule).mean())
        r['any'] = float((v != 'PASS').mean())
        rows.append(r)
    # the cohort's own pre-QC genotypes (03_annotate_af_norm_gt), the rules applied here, summed per sample
    keep = st[st.in_full_mainland == '1']
    group = {r.sample_id: ('control' if r.status == 'control' else ('deep' if r.deep == '1' else r.platform))
             for r in keep.itertuples(index=False)}
    rules = ('LowGQ', 'LowDP', 'ABB_outlier', 'ABB_NaN', 'any')
    sums, n_seen, names = None, 0, None

    def num(col):
        return pd.to_numeric(pd.Series(col).replace({'.': np.nan, '': np.nan, 'nan': np.nan, 'NaN': np.nan}),
                             errors='coerce').values

    for ch in sorted(set(fs.chrom), key=lambda c: int(c[3:])):
        vcf = a.cohort_vcf.replace('@CHROM@', ch)
        reg = Path(f'{ch}.regions.tsv')
        fs[fs.chrom == ch][['chrom', 'pos']].drop_duplicates().to_csv(reg, sep='\t', index=False, header=False)
        if names is None:
            allnames = subprocess.run([a.bcftools, 'query', '-l', vcf], check=True, capture_output=True,
                                      text=True).stdout.split()
            cols = [i for i, s in enumerate(allnames) if s in group]
            names = [allnames[i] for i in cols]
            sums = np.zeros((len(names), len(rules)))
        q = subprocess.Popen([a.bcftools, 'query', '-R', str(reg), '-f', '%ID[\t%GT:%GQ:%DP:%AF]\n', vcf],
                             stdout=subprocess.PIPE, text=True)
        for ln in q.stdout:
            f = ln.rstrip('\n').split('\t')
            if f[0] not in ids:
                continue
            x = np.array(f[1:], dtype=object)[cols]
            parts = np.array([v.split(':') for v in x])
            gt = parts[:, 0]
            alt = np.array([-1 if '.' in g else g.replace('|', '/').split('/').count('1') for g in gt])
            rr = rule_reasons(alt, num(parts[:, 1]), num(parts[:, 2]), num(parts[:, 3]))
            with np.errstate(invalid='ignore'):
                flags = np.vstack([np.nan_to_num(rr[k]).astype(bool) for k in rules[:4]])
            sums[:, :4] += flags.T
            sums[:, 4] += flags.any(axis=0)
            n_seen += 1
        if q.wait():
            sys.exit(f'ABORT: bcftools query failed on {vcf}')
    coh = pd.DataFrame(sums / max(n_seen, 1), columns=rules)
    coh.insert(0, 'n', n_seen)
    coh.insert(0, 'sample_id', names)
    coh.insert(0, 'arm', 'original')
    coh.insert(0, 'group', [group[s] for s in names])
    coh.insert(0, 'source', 'cohort')
    rec = pd.DataFrame(rows)
    per = pd.concat([rec, coh], ignore_index=True)
    per.to_csv(a.out_samples, sep='\t', index=False)
    summ = per.groupby(['source', 'group', 'arm'])[['LowGQ', 'LowDP', 'ABB_outlier', 'ABB_NaN', 'any']].median()
    summ['samples'] = per.groupby(['source', 'group', 'arm']).size()
    summ.reset_index().to_csv(a.out_summary, sep='\t', index=False)
    print(summ.round(4).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('variants')
    p.add_argument('--vcf', required=True)
    p.add_argument('--original', nargs='+', required=True, help='original fileset prefixes (full mainland first)')
    p.add_argument('--fixed', required=True)
    p.add_argument('--plan', required=True)
    for k in ('--out-long', '--out-summary', '--out-samples', '--out-replicates'):
        p.add_argument(k, required=True)
    p = sub.add_parser('reasons')
    p.add_argument('--vcf', required=True)
    p.add_argument('--fixed', required=True)
    p.add_argument('--plan', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--cohort-vcf', required=True, help='pre-QC cohort VCF with @CHROM@')
    p.add_argument('--bcftools', required=True)
    p.add_argument('--out-samples', required=True)
    p.add_argument('--out-summary', required=True)
    a = ap.parse_args()
    {'variants': cmd_variants, 'reasons': cmd_reasons}[a.cmd](a)


if __name__ == '__main__':
    main()
