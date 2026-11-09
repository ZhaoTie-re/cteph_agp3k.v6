#!/usr/bin/env python3
"""Stage 00 — one row per genome of the three cohorts: status, platform, cohort membership,
model sample sets, the CRAM and its index, and where an original GVCF survives.

Membership is read from the files the analyses read (each cohort's pheno.tsv for the random
model, the fixed-model .fam for the unrelated samples, tuning.rv's call-set .fam for the gene
tests), never re-derived. The platform label is the info sheet's WGS_Platform, verbatim.

CRAMs: `Cram_Path` is followed to the real file (the index sits beside it, not beside the
link). Links whose target directory is gone fall back to the same file name under
Control_BBJ with the index kept by hla.typing.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402


def cohort_sets():
    sets = {}
    for c in K.COHORTS:
        ph = pd.read_csv(K.MODEL_INPUTS / c / 'phenotype' / 'pheno.tsv', sep='\t', dtype=str)
        fixed = pd.read_csv(K.MODEL_INPUTS / c / 'genotype' / 'fixed_model' / f'{K.FIXED_PREFIX}.fam',
                            sep=r'\s+', header=None, dtype=str)[1]
        rv = pd.read_csv(K.TUNING_RV / c / '02.callset_filter' / 'filtered.fam', sep=r'\s+', header=None,
                         dtype=str)[1]
        sets[c] = {'pheno': dict(zip(ph.IID, ph.PHENO1)), 'fixed': set(fixed), 'rv': set(rv)}
    return sets


def resolve_cram(sid, path):
    """(cram, crai, source) or raises."""
    real = os.path.realpath(path)
    if os.path.exists(real) and os.path.exists(real + '.crai'):
        return real, real + '.crai', 'sheet'
    alt = K.CRAM_FALLBACK_DIR / f'{sid}.cram'
    crai = K.CRAI_FALLBACK_DIR / f'{sid}.cram.crai'
    if alt.exists() and crai.exists():
        return str(alt), str(crai), 'fallback'
    raise RuntimeError(f'{sid}: no readable CRAM + index for {path}')


def gvcf_dirs(ids):
    """sample -> the first GVCF folder (root order) holding all 22 autosomal GVCFs and indexes."""
    found, n_cand = {}, {}
    want = set(ids)
    for root in K.GVCF_ROOTS:
        for pat in ('*/gvcf.chr1-22', '*/*/gvcf.chr1-22'):
            for d in root.glob(pat):
                sid = d.parent.name
                if sid not in want:
                    continue
                n_cand[sid] = n_cand.get(sid, 0) + 1
                if sid in found:
                    continue
                if all((d / f'{sid}.chr{c}.g.vcf.gz').exists() and (d / f'{sid}.chr{c}.g.vcf.gz.tbi').exists()
                       for c in range(1, 23)):
                    found[sid] = str(d)
    return found, n_cand


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    info = pd.read_excel(K.INFO_XLSX, dtype=str).fillna('')
    sets = cohort_sets()
    full_ids = list(sets['full_mainland']['pheno'])
    for c in K.COHORTS:                      # nested cohorts: every narrower set inside full
        extra = set(sets[c]['pheno']) - set(full_ids)
        if extra:
            sys.exit(f'ABORT: {len(extra)} {c} samples are not in full_mainland')
    info = info[info.ID_JHRPv6.isin(full_ids)].drop_duplicates('ID_JHRPv6')
    if len(info) != len(full_ids):
        sys.exit(f'ABORT: info sheet holds {len(info)} of {len(full_ids)} full-mainland samples')
    unknown = sorted(set(info.WGS_Platform) - set(K.PLATFORM_SLUG))
    if unknown:
        sys.exit(f'ABORT: unknown WGS_Platform label(s) {unknown}')

    v6 = {ln.split('\t')[0] for f in K.V6_BATCH_LISTS for ln in f.read_text().splitlines() if ln.strip()}
    gv, n_cand = gvcf_dirs(full_ids)

    rows = []
    for r in info.itertuples(index=False):
        sid = r.ID_JHRPv6
        status = 'case' if r.Outcome == 'PH' else 'control' if r.Outcome == 'AGP3K' else ''
        if not status:
            sys.exit(f'ABORT: {sid} has Outcome {r.Outcome!r}')
        cram, crai, src = resolve_cram(sid, r.Cram_Path)
        row = {'sample_id': sid, 'status': status, 'platform': r.WGS_Platform,
               'platform_slug': K.PLATFORM_SLUG[r.WGS_Platform], 'target_depth': r.Target_Depth,
               'observed_depth': r.Observed_Depth, 'project': r.Project, 'sex': r.Sex,
               'deep': int(status == 'case' and r.WGS_Platform in K.DEEP_PLATFORMS),
               'cram': cram, 'crai': crai, 'cram_source': src,
               'batch': 'v6' if sid in v6 else ('v5' if (K.GVCF_ROOTS[0] / sid).exists() else 'v4_or_earlier'),
               'orig_gvcf_dir': gv.get(sid, ''), 'n_gvcf_candidates': n_cand.get(sid, 0)}
        for c in K.COHORTS:
            ph = sets[c]['pheno'].get(sid)
            row[f'in_{c}'] = int(ph is not None)
            if ph is not None and str(int(ph == '1')) != ('1' if status == 'case' else '0'):
                sys.exit(f'ABORT: {sid} is {status} in the info sheet but PHENO1={ph} in {c}')
            row[f'fixed_{c}'] = int(sid in sets[c]['fixed'])
            row[f'rv_{c}'] = int(sid in sets[c]['rv'])
        rows.append(row)
    df = pd.DataFrame(rows).sort_values(['status', 'platform', 'sample_id'])
    K.write_tsv(df, out / 'sample_table.tsv')

    # cohort x status x platform counts (gate G0's numbers)
    cnt = []
    for c in K.COHORTS:
        d = df[df[f'in_{c}'] == 1]
        for (st, pf), g in d.groupby(['status', 'platform']):
            cnt.append({'cohort': c, 'status': st, 'platform': pf, 'n': len(g),
                        'n_fixed_model': int(g[f'fixed_{c}'].sum()), 'n_rare_variant': int(g[f'rv_{c}'].sum())})
    cnt = pd.DataFrame(cnt)
    K.write_tsv(cnt, out / 'counts.tsv')
    deep = df[df.deep == 1]
    summary = {
        'n_samples': len(df), 'n_case': int((df.status == 'case').sum()),
        'n_control': int((df.status == 'control').sum()),
        'cram_fallback': int((df.cram_source == 'fallback').sum()),
        'deep_cases': len(deep),
        'deep_cases_by_platform': deep.platform.value_counts().to_dict(),
        'deep_cases_in_every_cohort': int(all(deep[f'in_{c}'].eq(1).all() for c in K.COHORTS)),
        'deep_cases_by_batch': deep.batch.value_counts().to_dict(),
        'deep_cases_with_original_gvcf': int((deep.orig_gvcf_dir != '').sum()),
        'controls_with_original_gvcf': int(((df.status == 'control') & (df.orig_gvcf_dir != '')).sum()),
    }
    (out / 'sample_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print(cnt.to_string(index=False))


if __name__ == '__main__':
    main()
