#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Collect every number the summary deck prints from the published tables of check/platform_depth
#           (read only) into stats.json. build_slides.js reads numbers only from stats.json, and deck_check.py checks
#           every number on the slides and in the notes against it.
# Inputs  : the results root and, optionally, a development copy of the report tables (environment, below);
#           references.json (literature, verified 2026-10-31); the design numbers of ../scripts/design_figures.py
#           (registry, results, run traces) and ../docs/design/figures.json
# Outputs : stats.json
# Environment (read by design_figures.py too, so the deck and docs/DESIGN.md use the same tables):
#   PD_RESULTS  the results root (default ../results)
#   PD_TABLES   optional: a flat folder of the report stage's tables in their final format (ds.*, array.*, controls.*,
#               followup.*, evidence.*), read instead of their published folders (e.g. ../.dev_report/tabs)
#   PD_FIGURES  optional: the component's figure families (default <PD_RESULTS>/figures); build_slides.js reads the
#               folder from stats.json
#   Without PD_TABLES the run is strict: every table must hold the arm plan's replicates and nothing else (the
#   superseded R = 10 arms must be gone), and the evidence tables may hold only the registry's checks.
# Usage   : /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 collect_stats.py
# ---------------------------------------------------------------------------
import glob
import json
import math
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
import design_figures as DF  # noqa: E402  (reads PD_RESULTS / PD_TABLES at import)

R = str(DF.RES)
FIGS = os.path.abspath(os.environ.get('PD_FIGURES') or os.path.join(R, 'figures'))
STRICT = DF.TABLES is None
CK = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
PLATS = ['DNBSeq-T7 30x', 'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x']
SLUG = {'t7': 'DNBSeq-T7 30x', 'novaseq': 'NovaSeq 30x', 'g400rs30': 'DNBSeq-G400RS 30x',
        'g400rs15': 'DNBseq-G400RS 15x'}
SOURCES = ['nagahama', 'bbj', 'htlv1', 'acc', 'lungca']
WARN = []


def tsv(*p, **kw):
    return pd.read_csv(os.path.join(R, *p), sep='\t', **kw)


def js(*p):
    return json.load(open(os.path.join(R, *p)))


def rtsv(name, **kw):
    return pd.read_csv(DF.report_file(name), sep='\t', **kw)


def rjs(name):
    return json.load(open(DF.report_file(name)))


def clean(x):
    """JSON-safe: NaN -> None, numpy -> python."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if hasattr(x, 'item'):
        x = x.item()
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def check(ok, msg):
    """a consistency rule of the final results: an error in the strict run, a warning in a development run"""
    if not ok:
        if STRICT:
            sys.exit(f'collect_stats: {msg}')
        WARN.append(msg)


S = {}
# ── the design numbers first: scripts/design_figures.py numbers() (registry, results, run traces) ───────────────
N = DF.numbers()
REPS = N['arms']['replicate_names']                    # the arm plan's replicates (r01 ...)
ARMS = ['original', 'baseline'] + REPS
S['roots'] = {'results': R, 'tables': str(DF.TABLES) if DF.TABLES else None, 'figures': FIGS,
                'mode': 'final' if STRICT else 'development'}
S['imputation'] = N.pop('imputation')                  # the one key every slide takes the imputation tool from
S['design'] = N

# ── samples ─────────────────────────────────────────────────────────────────
cnt = tsv('00.samples', 'counts.tsv')
S['samples'] = {c: {r.platform: {'status': r.status, 'n': r.n, 'n_fixed': r.n_fixed_model, 'n_rv': r.n_rare_variant}
                    for r in cnt[cnt.cohort == c].itertuples()} for c in CK}
S['n_case'] = {c: int(cnt[(cnt.cohort == c) & (cnt.status == 'case')].n.sum()) for c in CK}
S['n_control'] = {c: int(cnt[(cnt.cohort == c) & (cnt.status == 'control')].n.sum()) for c in CK}
ss = js('00.samples', 'sample_summary.json')
S['deep'] = {'n': ss['deep_cases'], 'by_platform': ss['deep_cases_by_platform'], 'cram_fallback': ss['cram_fallback']}
S['t7_share_of_cases'] = {c: S['samples'][c]['DNBSeq-T7 30x']['n'] / S['n_case'][c] for c in CK}
S['g400rs_cases'] = S['samples']['full_mainland']['DNBSeq-G400RS 30x']['n'] + \
    S['samples']['full_mainland']['DNBseq-G400RS 15x']['n']

# ── depth ───────────────────────────────────────────────────────────────────
pdp = tsv('01.depth', 'platform_depth.tsv')
S['depth'] = {r.platform: {'n': r.n, 'readlen': r.readlen_median, 'depth': r.depth_median, 'q25': r.depth_q25,
                           'q75': r.depth_q75, 'dup': r.dup_share_median}
              for r in pdp[pdp.cohort == 'full_mainland'].itertuples()}
dsm = js('01.depth', 'depth_summary.json')
S['target'] = {'depth': dsm['target_depth'], 'fraction_median': dsm['fraction_median'],
               'fraction_by_platform': dsm['fraction_by_platform'],
               'local_ratio': dsm['local_ratio_deep_after_over_controls']}
loc = tsv('01.depth', 'local_depth.tsv')
t7f = N['downsampled']['t7_fraction']
lr = loc.median_t7 * t7f / loc.median_hiseqx
S['target']['local_ratio_t7'] = [float(lr.min()), float(lr.median()), float(lr.max())]
S['target']['local_ratio_t7_before'] = [float((loc.median_t7 / loc.median_hiseqx).min()),
                                        float((loc.median_t7 / loc.median_hiseqx).max())]

# ── gates ───────────────────────────────────────────────────────────────────
G = {k: N['gates'][k] for k in ('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G8', 'G9', 'G10')}
g3 = js('_gates', 'G3.json')
G['G3']['pass'] = g3['pass']
G['G6_pass'] = js('_gates', 'downsampling', 'G6.json')['pass']
G['G7'] = {'swap_records': N['arms']['swap_records']}
# the cases swapped into each cohort's files (the same in every arm; the gene analyses hold their own samples)
sw = {}
for f in sorted(glob.glob(os.path.join(R, '04.genotypes', 'downsampling', '*', '*', 'swap.*.json'))):
    arm, coh, kind = f.split(os.sep)[-3], f.split(os.sep)[-2], os.path.basename(f).split('.')[1]
    d = json.load(open(f))
    sw.setdefault(coh, {}).setdefault(kind, set()).add(int(d['samples_swapped']))
    check(d.get('untouched_identical') is True, f'{f}: an untouched genome changed')
check(all(len(v) == 1 for k in sw.values() for v in k.values()), f'swapped cases differ between arms: {sw}')
G['G7']['swapped'] = {c: {k: v.pop() for k, v in sw[c].items()} for c in CK}
G['pass'] = N['gates']['pass']
ap = tsv('03.calls', 'downsampling', 'arm_plan.tsv')
fid = ap[ap.role == 'fidelity']
G['fidelity'] = {'n': int(fid.sample_id.nunique()), 'platforms': sorted(fid.platform.unique().tolist())}
S['gates'] = G
# G6 in one line (the appendix gates table): comparisons, the lowest concordance, the largest call-rate difference
S['g6_summary'] = {'comparisons': len(G['G6']), 'conc_min': min(v['concordance'] for v in G['G6'].values()),
                   'callrate_maxdiff': max(abs(v['callrate_ref'] - v['callrate_test']) for v in G['G6'].values()),
                   'controls': int(G['fidelity']['n'])}

# ── replicates: the user's decision (three for every platform); the pilot's measured costs ──────────────────
P = N['pilot']
S['pilot'] = {'R': P['R'], 'R_rule': P['R_rule'], 'rule': P['rule'], 'projection': P['projection'],
              'core_h_per_task': P['core_h_per_task'], 'tasks_measured': P['tasks_measured']}
m = re.search(r'R in \[(\d+), (\d+)\] with cost <= ([\d.]+) core-h and wall <= ([\d.]+) h at (\d+) cores', P['rule'])
S['pilot'].update({'R_min': int(m.group(1)), 'R_max': int(m.group(2)), 'cost_limit': float(m.group(3)),
                   'wall_limit': float(m.group(4)), 'cores': int(m.group(5))})
S['pilot']['projection_R'] = next(x for x in P['projection'] if x['R'] == P['R'])
S['pilot']['projection_rule'] = next(x for x in P['projection'] if x['R'] == P['R_rule'])
# the measured cost per down-sampled genome and arm: the per-task cost of the per-genome processes of an arm
per_arm = ['EXTRACT_READS', 'SUBSAMPLE', 'CALL_GVCF', 'HLA_SUBSAMPLE', 'HLAHD']
S['pilot']['core_h_per_genome_arm'] = sum(P['core_h_per_task'][k] for k in per_arm if k in P['core_h_per_task'])
S['pilot']['hlahd_core_h_per_task'] = P['core_h_per_task']['HLAHD']

# ── down-sampling: genotypes (genome-level, the same in every cohort) ─────────────────────────────────────────
DS = os.path.join('07.downsampling', 'downsampling')
rep = tsv(DS, 'concordance', 'concordance.replicates.tsv')
check(set(rep.replicates.dropna().astype(int)) <= {len(REPS)},
      f'concordance.replicates.tsv summarizes {sorted(set(rep.replicates.dropna().astype(int)))} replicates, the arm plan {len(REPS)}')
S['ds_geno'] = {}
for plat in ['all'] + sorted(rep.platform.unique()):
    a = rep[(rep.comparison == 'replicate_vs_baseline') & (rep.role == 'deep') & (rep.platform == plat) & (rep.kind == 'all')]
    if not len(a):
        continue
    a = a.iloc[0]
    S['ds_geno'][plat] = {'conc_median': a.concordance_median, 'conc_min': a.concordance_min,
                          'conc_max': a.concordance_max, 'nonref_min': a.nonref_concordance_min,
                          'het_median': a.het_retained_median, 'het_min': a.het_retained_min,
                          'het_max': a.het_retained_max, 'callrate_base': a.callrate_ref_median,
                          'callrate_min': a.callrate_test_min, 'callrate_max': a.callrate_test_max,
                          'callrate_median': a.callrate_test_median}
S['ds_geno_platforms'] = [p for p in S['ds_geno'] if p != 'all']
# QC-rule failures per cohort (design numbers: median share over each group's genomes)
S['qc_fail'] = N['qc_fail']
qr = rtsv('ds.qc_reasons.tsv')
extra = sorted(set(qr.arm) - set(ARMS) - {'origgvcf'})
check(not extra, f'ds.qc_reasons.tsv holds arms outside the arm plan: {extra}')
qq = qr[qr.arm.isin(REPS) & (qr.source == 'recalled') & (qr.group == 'deep') & (qr.platform != 'all')]
S['qc_rule_reps'] = {c: {p: {k: float(qq[(qq.cohort == c) & (qq.platform == p)][k].median())
                             for k in ('LowGQ', 'LowDP', 'ABB_outlier', 'any')}
                         for p in sorted(qq[qq.cohort == c].platform.unique())} for c in CK}
co = qr[(qr.source == 'cohort') & (qr.group == 'control') & (qr.arm == 'original')]
S['qc_rule_controls'] = {r.cohort: {k: float(getattr(r, k)) for k in ('LowGQ', 'LowDP', 'ABB_outlier', 'any')}
                         for r in co.itertuples()}
# the HiSeq X fidelity controls re-processed at f = 1 (their own share, not the controls' median)
fq = qr[(qr.source == 'recalled') & (qr.group == 'fidelity') & (qr.platform == 'all') & (qr.arm == 'baseline')]
S['qc_fail_fidelity'] = {r.cohort: float(r.any) for r in fq.itertuples()}
# how many of the fidelity controls each cohort holds (the cohorts are nested; full mainland holds every one)
S['qc_fail_fidelity_n'] = {r.cohort: int(r.samples) for r in fq.itertuples()}
check(S['qc_fail_fidelity_n'].get('full_mainland') == G['fidelity']['n'],
      f"ds.qc_reasons.tsv: {S['qc_fail_fidelity_n']} fidelity controls in full mainland, the arm plan {G['fidelity']['n']}")
# the down-sampled cases whose published GVCFs were re-genotyped (gate G6), by platform (full mainland)
og = qr[(qr.cohort == 'full_mainland') & (qr.source == 'recalled') & (qr.group == 'deep') & (qr.arm == 'origgvcf')
        & (qr.platform != 'all')]
S['orig_gvcf_cases'] = {r.platform: int(r.samples) for r in og.itertuples()}
hs = tsv(DS, 'hla', 'hla.summary.tsv')
hr = hs[hs.comparison == 'replicate_vs_baseline']
if 'arm' in hr.columns:
    check(set(hr.arm) <= set(REPS), f'hla.summary.tsv holds arms {sorted(set(hr.arm) - set(REPS))} outside the arm plan')
    hr = hr[hr.arm.isin(REPS)]
S['ds_hla'] = {g: {'conc_min': hr[hr.gene == g].allele_concordance.min(),
                   'conc_max': hr[hr.gene == g].allele_concordance.max(),
                   'reads_base': hr[hr.gene == g].reads_ref_median.iloc[0],
                   'reads_rep_min': hr[hr.gene == g].reads_test_median.min(),
                   'reads_rep_max': hr[hr.gene == g].reads_test_median.max()} for g in hr.gene.unique()}
S['ds_hla_all'] = {'conc_min': hr.allele_concordance.min(), 'conc_max': hr.allele_concordance.max(),
                   'genes': int(hr.gene.nunique())}

# ── the signals: evidence matrix ───────────────────────────────────────────
ev = rtsv('evidence.signals.tsv')
lg = rtsv('evidence.long.tsv')
loo = rtsv('evidence.loo.tsv')
# the deck leaves the HLA residue position (DQB1 position 57) out (user 2026-11-05: redundant with DQB1:57:V); the
# pipeline still tests and reports it. Every count on the slides is over the remaining signals.
POS_SIGS = sorted(set(ev[ev.kind == 'hla_position'].signal))
lg_pos = lg[lg.signal.isin(POS_SIGS)]
ev, lg, loo = ev[~ev.signal.isin(POS_SIGS)], lg[~lg.signal.isin(POS_SIGS)], loo[~loo.signal.isin(POS_SIGS)]
rules = rtsv('evidence.rules.tsv')
extra = sorted((set(lg.check) | set(rules.check)) - set(DF.CK.CHECK))
check(not extra, f'the evidence tables hold checks outside the registry (scripts/checks.py): {extra}')
S['rules'] = dict(zip(rules.check, rules.rule))
# the residue position is not in the deck: its threshold drops out of the discovery rule as the slides quote it
check('; positions Bonferroni' in S['rules']['discovery'], 'the discovery rule no longer names the positions')
S['rules']['discovery'] = S['rules']['discovery'].replace('; positions Bonferroni', '')
S['signals'] = [clean(r._asdict()) for r in ev.itertuples(index=False)]
S['long'] = {}
for r in lg.itertuples(index=False):
    S['long'].setdefault(r.signal, {}).setdefault(r.check, {})[r.cohort] = {
        'reading': r.reading, 'metric': None if pd.isna(r.metric) else r.metric,
        'value': None if pd.isna(r.value) else r.value}
# the readings of each check (a reading = a signal and cohort where the check applies), and how many are consistent
app = lg[lg.reading != 'not applicable']
S['check_readings'] = {k: {'readings': int(len(v)), 'consistent': int((v.reading == 'consistent').sum()),
                           'flag': int((v.reading == 'flag').sum())} for k, v in app.groupby('check')}
S['loo'] = {}
for r in loo.itertuples(index=False):
    S['loo'].setdefault(r.signal, {}).setdefault(r.check, {}).setdefault(r.cohort, {})[r.left_out] = {
        'z': r.z_observed, 'z_exp': r.z_expected}
chk = lg[lg.check != 'discovery']
S['matrix'] = {'cells': {k: int(v) for k, v in chk.reading.value_counts().items()}, 'readings': int(len(chk)),
               'discovery_readings': int((lg.check == 'discovery').sum())}
S['loo_points'] = int(len(loo))
S['flag_cells'] = int((chk.reading == 'flag').sum())
S['flag_signals'] = sorted(ev[ev['flags'] != 'none'].signal.tolist())
S['n_signals'] = int(len(ev))
S['checks'] = [c for c in ev.columns if c not in ('signal', 'kind', 'discovery', 'flags')]
S['n_checks'] = len(S['checks'])
# signals contradicted by the checks that can contradict (depth, platform heterogeneity, arrays, control sources)
contra = ['depth', 'platform_Q', 'platform_RE', 'array', 'sources_Q', 'sources_loo']
S['n_contradicted'] = int(chk[chk.check.isin(contra) & (chk.reading == 'flag')].signal.nunique())

# discovery P per signal and cohort (the ORIGINAL arm = the published result)
disc = lg[lg.check == 'discovery']
S['discovery'] = {s: {r.cohort: {'reading': r.reading, 'p': r.value, 'metric': r.metric} for r in g.itertuples()}
                  for s, g in disc.groupby('signal')}
gw = ev[ev.kind == 'gwas'].signal.tolist()
S['gw_loci'] = [g for g in gw if any(S['discovery'][g][c]['reading'] == 'genome-wide' for c in CK)]
S['sg_only'] = [g for g in gw if all(S['discovery'][g][c]['reading'] == 'suggestive' for c in CK)]
sgp = [S['discovery'][g][c]['p'] for g in S['sg_only'] for c in CK]
S['sg_p_range'] = [min(sgp), max(sgp)]

# down-sampling per signal (GWAS leads, HLA): |z| original and replicate median from the metric text
dep = lg[lg.check == 'depth']
S['ds_signal'] = {}
for r in dep.itertuples():
    m1 = re.search(r'\|z\| original ([\d.]+), replicate median ([\d.]+)', r.metric or '')
    m2 = re.search(r'P < 1e-5: (\d+)/(\d+), < 5e-8: (\d+)/(\d+)', r.metric or '')
    m3 = re.search(r'significant(?: at GEC)?: (\d+)/(\d+)', r.metric or '')
    m4 = re.search(r'case carriers (\d+) -> (\d+)-(\d+)', r.metric or '')
    d = {'metric': r.metric, 'dz': r.value}
    if m1:
        d.update(z_orig=float(m1.group(1)), z_med=float(m1.group(2)))
    if m2:
        d.update(reps_sg=int(m2.group(1)), reps=int(m2.group(2)), reps_gw=int(m2.group(3)))
    if m3:
        d.update(reps_sig=int(m3.group(1)), reps=int(m3.group(2)))
    if m4:
        d.update(carriers=int(m4.group(1)), carriers_min=int(m4.group(2)), carriers_max=int(m4.group(3)))
    if 'reps' in d:
        check(d['reps'] == len(REPS), f'{r.signal} {r.cohort}: the depth reading counts {d["reps"]} replicates')
    S['ds_signal'].setdefault(r.signal, {})[r.cohort] = d
drops = {c: [S['ds_signal'][s][c]['dz'] for s in gw] for c in CK}   # unrounded (reading value)
S['ds_gwas'] = {'max_drop': round(max(max(v) for v in drops.values()), 2),
                'max_drop_by_cohort': {c: round(max(v), 2) for c, v in drops.items()},
                'leads': len(gw), 'flags': {c: int(((dep.cohort == c) & (dep.reading == 'flag')).sum()) for c in CK}}

# DQB1:57:V across the replicates
v = rtsv('ds.hla_dqb1_57v.tsv')
S['ds_v57'] = {}
for r in v.to_dict('records'):
    ors = [r[f'or_{a}'] for a in REPS]
    S['ds_v57'][r['cohort']] = {'or_orig': r['or_original'], 'or_min': min(ors), 'or_max': max(ors),
                                'p_orig': r['p_original'], 'reps_gec': r['reps_gec'], 'reps': r['reps'],
                                'pg_orig': r['p_over_gec_original'], 'pg_min': r['p_over_gec_rep_min'],
                                'pg_max': r['p_over_gec_rep_max'], 'pg_base': r['p_over_gec_baseline'],
                                'pg': {a: r[f'p_over_gec_{a}'] for a in REPS},
                                'flip_arms': [a for a in REPS if r[f'p_over_gec_{a}'] >= 1]}
S['ds_n3'] = rtsv('ds.hla_n3.tsv').to_dict('records')
dsj = rjs('ds.summary.json')
check(dsj.get('replicates') == REPS, f'ds.summary.json replicates {dsj.get("replicates")} differ from the arm plan {REPS}')
S['ds_v57_copies'] = dsj['dqb1_57v_copies']
fv = rtsv('followup.hla_v57.tsv')
extra = sorted(set(fv.arm) - set(ARMS))
check(not extra, f'followup.hla_v57.tsv holds arms outside the arm plan: {extra}')
fv = fv[fv.arm.isin(ARMS) & (fv.cohort == 'full_mainland')]
# the follow-up's own criterion (signal_followup.hla_v57), within the arm plan: the 57 V dosage differs between two
# arms, a missing dosage counted as a value (make_figures.v57_trace draws the same cases)
S['v57_traced'] = sorted(s for s, x in fv.groupby('sample_id') if x.v57_dosage.astype(str).nunique() > 1)
fv = fv[fv.sample_id.isin(S['v57_traced'])]
S['v57_trace'] = clean(fv.to_dict('records'))
# the platform of each traced case: the deck names a case by its platform, never by its sample ID
_plat = tsv('00.samples', 'sample_table.tsv', dtype=str).set_index('sample_id').platform
S['v57_traced_platform'] = {s: str(_plat.get(s, '')) for s in S['v57_traced']}
# the down-sampled cases whose 57 V dosage changes in a replicate (or whose call turns ambiguous)
chg = fv[fv.arm.isin(REPS)]
S['v57_changed'] = sorted(chg[(chg.change_vs_baseline.fillna(0) != 0)].sample_id.unique().tolist())
S['v57_ambiguous'] = sorted(chg[chg.v57_dosage.isna()].sample_id.unique().tolist())
# how many down-sampled cases carry each baseline DQB1 allele of a traced case (the baseline calls of the re-typing)
hc = tsv(DS, 'hla', 'hla.calls_long.tsv', dtype=str)
hb = hc[(hc.gene == 'DQB1') & (hc.comparison == 'replicate_vs_baseline')].drop_duplicates('sample_id')
held = hb.alleles_ref.fillna('').str.split('/').apply(lambda x: {a for a in x if a})
S['v57_baseline_typed'] = int(len(hb))
# every DQB1 call of a replicate that differs from the baseline call (allele sets), and the cases behind them
rc = hc[(hc.gene == 'DQB1') & (hc.comparison == 'replicate_vs_baseline') & hc.arm.isin(REPS)]
S['dqb1_changed'] = {'calls': int((rc.identical != '1').sum()), 'of_calls': int(len(rc)),
                     'cases': int(rc[rc.identical != '1'].sample_id.nunique())}
S['v57_allele_carriers'] = {}
for s in S['v57_traced']:
    h1 = held[hb.sample_id == s]
    for a in (sorted(h1.iloc[0]) if len(h1) else []):
        S['v57_allele_carriers'][a] = int(held.apply(lambda x, a=a: a in x).sum())

# Tier 1-2 gene sets across the replicates (CMC and SKAT-O)
dg = rtsv('ds.genes.tsv')
S['ds_genes'] = {}
for r in dg.to_dict('records'):
    S['ds_genes'][f"{r['gene']}|{r['stratum']}|{r['method']}|{r['cohort']}"] = {
        'p_orig': r['pvalue_original'], 'p_min': min(r[f'pvalue_{a}'] for a in REPS),
        'p_max': max(r[f'pvalue_{a}'] for a in REPS),
        'case_orig': r['carriers_case_original'], 'case_min': r['carriers_case_rep_min'],
        'case_max': r['carriers_case_rep_max'], 'ctrl_orig': r['carriers_control_original'],
        'ctrl_min': min(r[f'carriers_control_{a}'] for a in REPS),
        'ctrl_max': max(r[f'carriers_control_{a}'] for a in REPS), 'reps_bonf': r['reps_bonferroni'], 'reps': r['reps'],
        # the replicates (1-based) whose case carriers differ from the original
        'reps_case_changed': [i + 1 for i, a in enumerate(REPS)
                              if r[f'carriers_case_{a}'] != r['carriers_case_original']]}
S['ds_genes_changed'] = sorted({k.split('|')[0] for k, d in S['ds_genes'].items()
                                if (d['case_min'], d['case_max']) != (d['case_orig'], d['case_orig'])
                                or (d['ctrl_min'], d['ctrl_max']) != (d['ctrl_orig'], d['ctrl_orig'])})

# exposure: the share of the case ALT copies the down-sampled cases carry at each lead, per cohort
exp = rtsv('ds.exposure.tsv')
psig = tsv('10.platform', '_summary', 'platform.signals.tsv')
leads = psig[(psig.cohort == 'full_mainland') & (psig.role.isin(['Lead', 'HLA residue']))][['signal', 'marker']]
S['lead_marker'] = dict(zip(leads.signal, leads.marker))
ex = exp[exp.variant_id.isin(leads.marker)]
# ... and the median share over every GWAS signal variant (leads, alternative leads, credible-set members)
eg = exp[exp.kind == 'gwas']
S['exposure_median'] = {c: float(eg[eg.cohort == c].downsampled_share_of_case_alt_copies.median()) for c in CK}
S['exposure_variants'] = int(eg[eg.cohort == 'full_mainland'].variant_id.nunique())
S['exposure'] = {}
for c in CK:
    e = ex[ex.cohort == c]
    tot = e.case_alt_copies.where(e.case_alt_copies > 0)
    S['exposure'][c] = {k: [float((e[f'{k}_alt_copies'] / tot).min()), float((e[f'{k}_alt_copies'] / tot).max())]
                        for k in ('downsampled', 'deep', 't7')}

# ── platform strata (unchanged by the down-sampling) ────────────────────────
S['platform'] = {}
for r in psig[psig.role.isin(['Lead', 'HLA residue'])].itertuples():
    S['platform'].setdefault(r.signal, {})[r.cohort] = {
        'or': r.or_pooled, 'p_Q': r.p_Q_empirical, 'I2': r.I2, 'or_re': r.or_re, 're_l95': r.or_re_l95,
        're_u95': r.or_re_u95, 'p_case_only': r.p_case_only,
        **{f'or_{k}': getattr(r, f'or_only_{k}') for k in SLUG},
        **{f'ci_{k}': getattr(r, f'or_only_{k}_ci') for k in SLUG},
        **{f'orw_{k}': getattr(r, f'or_without_{k}') for k in SLUG},
        **{f'pw_{k}': getattr(r, f'p_without_{k}') for k in SLUG},
        **{f'pdiff_{k}': getattr(r, f'p_diff_{k}') for k in SLUG}}
S['platform_null'] = {c: js('10.platform', c, 'platform.summary.json') for c in CK}
m = re.search(r'([\d.]+) \(threshold ([\d.]+)\)', S['long']['FTO']['platform_Q']['full_mainland']['metric'])
S['q_threshold'] = float(m.group(2))
qs = [g for g in gw + ['DQB1:57:V'] if S['long'][g]['platform_Q']['full_mainland']['reading'] != 'not applicable']
S['q_markers'] = qs
S['q_min'] = {'platform': min(S['platform'][g][c]['p_Q'] for g in qs for c in CK)}
# per-stratum allele frequencies and ORs of the two flagged leads (SAIGE allele frequencies), every cohort
S['strata_af'] = {}
for c in CK:
    st = tsv('10.platform', c, 'platform.strata.tsv')
    for sig in ('ADARB2', 'LINC01146'):
        mk = S['lead_marker'][sig]
        for r in st[(st.variant_id == mk)].to_dict('records'):
            S['strata_af'].setdefault(sig, {}).setdefault(c, {})[r['subset']] = {
                k: r[k] for k in ('af_case', 'af_ctrl', 'n_case', 'or', 'or_l95', 'or_u95', 'p')}

# ── arrays, every cohort ────────────────────────────────────────────────────
asum = rjs('array.summary.json')
S['array'] = {c: N['array'][c] for c in CK}
S['array_gate'] = N['G10']
# DQB1:57:V on the arrays: HIBAG with its published Asian model (array_hibag.py; user 2026-11-05)
S['array_hibag'] = {c: js('11.array', 'cohorts', c, '07.hibag', 'summary.json') for c in CK}
# the two array data sets Q3 reads (user 2026-11-05: the Q3 schematic): the QC'd arrays and the Minimac4-imputed arrays
_aq = DF.ROOT.parents[1].parent / 'cteph_agp3k.v5' / 'array' / 'results' / '10.variant_qc_filtered' / 'cteph_agp3k.array.sqc.vqc'
_fam = pd.read_csv(f'{_aq}.fam', sep=r'\s+', header=None)
_ai = DF.ROOT.parents[1].parent / 'cteph_agp3k' / 'array_imputed_rev1' / 'results' / '04.sample_qc'
_log = (_ai / 'cteph_agp3k.imputed_array.sqc.log').read_text()
S['array_sets'] = {
    'qc': {'samples': int(len(_fam)), 'cases': int((_fam[5] == 2).sum()), 'controls': int((_fam[5] == 1).sum()),
           'variants': int(sum(1 for _ in open(f'{_aq}.bim')))},
    'imputed': {'samples': int(re.search(r'(\d+) samples .* remaining after main', _log).group(1)),
                'variants': int(re.search(r'(\d+) variants loaded', _log).group(1)), 'tool': 'Minimac v4.1.4',
                'panel': 'AGP3K whole-genome reference panel', 'chromosomes': '1–22'}}
check(S['array_sets']['qc']['samples'] == S['array_sets']['imputed']['samples'], 'QC\'d and imputed arrays differ in people')
# HIBAG and HLA-HD side by side (user 2026-11-05: a slide on the HIBAG method against HLA-HD): the DQB1 model's own
# facts (assets/hibag/DQB1_model.json, written from the pinned model by scripts/hibag_model_summary.R) and HLA-HD's
# reference and read filter as the hla.typing component records them
S['hibag_model'] = json.load(open(DF.ROOT / 'assets' / 'hibag' / 'DQB1_model.json'))
_rd = (DF.ROOT.parent.parent / 'hla.typing' / 'README.md').read_text()
_m = re.search(r'IPD-IMGT/HLA ([\d.]+), ([\d,]+) alleles', _rd)
S['hlahd_facts'] = {'imgt': _m.group(1), 'alleles': int(_m.group(2).replace(',', '')),
                    'min_read': int(re.search(r'`MinReadLength` \| (\d+)', _rd).group(1)), 'version': '1.7.1'}
# the numbers the HLA-HD slide of hla.typing.v2 quotes (its deck/stats.json): genomes typed and how each 2-field name
# was settled
_hv = json.load(open(DF.ROOT.parent.parent / 'hla.typing.v2' / 'deck' / 'stats.json'))
S['hlahd_basis'] = {**_hv['basis'], 'genomes': int(_hv['cohort']['n']), 'validation_ok': int(_hv['cohort']['validation_ok'])}
asig = rtsv('array.signals.tsv')
S['array_signal'] = {c: {r.signal: {'evidence': r.strongest_evidence, 'verdict': r.array_verdict,
                                    'recovered': None if pd.isna(r.case_carriers_recovered) else r.case_carriers_recovered}
                         for r in asig[asig.cohort == c].itertuples()} for c in CK}
aev = rtsv('array.evidence.tsv')
cols = ['cohort', 'signal', 'kind', 'marker', 'evidence', 'array_marker', 'info_r2', 'proxy_r2', 'or_wgs', 'p_wgs',
        'or_array', 'p_array', 'dz', 'z_wgs', 'z_array', 'verdict', 'or_allele', 'conc_cases', 'conc_controls', 'nonref_conc_cases',
        'nonref_conc_controls', 'carriers', 'carriers_wgs',
        'carriers_found', 'expected_recovery',
        'carriers_expected', 'mean_dosage', 'carriers_extra', 'controls_wgs', 'controls_found', 'posterior_median']
S['array_rows'] = clean(aev[[c for c in cols if c in aev.columns]].to_dict('records'))
# how every signal is matched (user 2026-11-05: one marker per signal: typed, else imputed with R2 >= 0.8, else an LD
# proxy; DQB1:57:V with HIBAG). The imputation R2 is one value per marker (the imputation ran once); the match is the
# same in the three cohorts (checked). Gene drivers (user 2026-11-05): full mainland only, the WGS case carriers looked
# up on the imputed arrays (no verdict)
S['array_match'] = []
_af = aev[aev.cohort == 'full_mainland']
for _r in asig[asig.cohort == 'full_mainland'].itertuples():
    _e = _af[_af.signal == _r.signal]
    _ck = CK if _r.kind != 'gene' else ['full_mainland']
    _by = {c: asig[(asig.cohort == c) & (asig.signal == _r.signal)].strongest_evidence.iloc[0] for c in _ck}
    check(len(set(_by.values())) == 1, f'{_r.signal}: the array match differs between cohorts')
    check(_r.kind != 'gene' or not len(asig[(asig.cohort != 'full_mainland') & (asig.signal == _r.signal)]),
          f'{_r.signal}: gene carriers outside full mainland')
    _row = {'signal': _r.signal, 'kind': _r.kind, 'evidence': _r.strongest_evidence,
            'verdict_by_cohort': {c: asig[(asig.cohort == c) & (asig.signal == _r.signal)].array_verdict.iloc[0]
                                  for c in _ck}}
    if _r.kind == 'gene':
        _num = lambda v, f=int: None if pd.isna(v) else f(v)  # noqa: E731
        _row['drivers'] = [{'marker': m.marker.split(' ', 1)[1], 'evidence': m.evidence,
                            'carriers_wgs': _num(m.carriers_wgs), 'carriers_found': _num(m.carriers_found),
                            'expected_recovery': _num(m.expected_recovery, float),
                            'carriers_expected': _num(m.carriers_expected, float), 'mean_dosage': _num(m.mean_dosage, float),
                            'carriers_extra': _num(m.carriers_extra), 'controls_wgs': _num(m.controls_wgs),
                            'info_r2': _num(m.info_r2, float), 'conc_cases': _num(m.conc_cases, float),
                            'conc_controls': _num(m.conc_controls, float), 'nonref_conc_cases': _num(m.nonref_conc_cases, float),
                            'controls_found': _num(m.controls_found), 'verdict': m.verdict} for m in _e.itertuples()]
    else:
        _x = _e.iloc[0]
        _row.update({'wgs_marker': str(_x.marker), 'array_marker': None if pd.isna(_x.array_marker) else str(_x.array_marker),
                     'info_r2': None if pd.isna(_x.get('info_r2', float('nan'))) else float(_x.info_r2),
                     'proxy_r2': None if pd.isna(_x.get('proxy_r2', float('nan'))) else float(_x.proxy_r2),
                     'conc_cases': None if pd.isna(_x.get('conc_cases', float('nan'))) else float(_x.conc_cases),
                     'conc_controls': None if pd.isna(_x.get('conc_controls', float('nan'))) else float(_x.conc_controls),
                     'nonref_conc_cases': None if pd.isna(_x.get('nonref_conc_cases', float('nan'))) else float(_x.nonref_conc_cases)})
        if _r.strongest_evidence == 'imputed (HIBAG)':
            _row['hibag'] = {'model_snps_on_array': int(S['array_hibag']['full_mainland']['run']['model_snps_on_array']),
                             'model_snps': int(S['array_hibag']['full_mainland']['run']['model_snps'])}
    S['array_match'].append(_row)
_g = [r for r in S['array_match'] if r['kind'] == 'gwas']
S['array_match_counts'] = {e: sum(r['evidence'] == e for r in _g) for e in ('typed', 'imputed', 'LD proxy', 'none')}
_dr = [d for r in S['array_match'] if r['kind'] == 'gene' for d in r['drivers']]
_dc = [d for d in _dr if d['verdict'] == 'counted']
_gc = asum['gene_carriers']
# recovery of WGS carriers on the imputed arrays (full mainland; array_cohorts.py step 06), over the rare untyped
# variants within 500 kb of the drivers (array_crosscheck's driver windows, not genome-wide):
# the controls are in the panel (near-complete recovery); the cases' recovery by number of control carriers
_rb = tsv('11.array', 'cohorts', 'full_mainland', '06.imputed', 'recovery_by_group.tsv')
_rbc = _rb[_rb.group == 'controls']
S['array_drivers'] = {'n': len(_dr), 'counted': len(_dc), 'genes': len({r['signal'] for r in S['array_match'] if r['kind'] == 'gene'}),
                      'case_wgs': sum(d['carriers_wgs'] for d in _dc), 'case_found': sum(d['carriers_found'] for d in _dc),
                      'case_expected': sum(d['carriers_expected'] for d in _dc),
                      'case_extra': sum(d['carriers_extra'] for d in _dc),
                      'recovery_range': [min(d['expected_recovery'] for d in _dc), max(d['expected_recovery'] for d in _dc)],
                      'not_in_imputed': [f"{r['signal']} {d['marker']}" for r in S['array_match'] if r['kind'] == 'gene'
                                         for d in r['drivers'] if d['verdict'] != 'counted'],
                      'controls_recovery_nearby': float(_rbc.recovered.sum() / _rbc.wgs_carriers.sum()),
                      'reference_variants': int(_rb[_rb.group == 'cases'].variants.sum()),
                      'window_kb': int(round(max((int(e) - int(b)) for _c, b, e, _g in
                                                 (l.split('\t') for l in open(DF.ROOT.parent / 'array_crosscheck' / 'results' /
                                                  '06.imputed' / 'drivers.range').read().split('\n') if l)) / 2 / 1000 / 100) * 100),
                      'case_recovery_by_bin': {r.bin: float(r.recovery) for r in _rb[_rb.group == 'cases'].itertuples()},
                      'reference_variants_by_bin': {r.bin: int(r.variants) for r in _rb[_rb.group == 'cases'].itertuples()},
                      'case_false_carrier_rate': [float(_rb[_rb.group == 'cases'].false_carrier_rate.min()),
                                                  float(_rb[_rb.group == 'cases'].false_carrier_rate.max())]}
check(S['array_drivers']['case_found'] == _gc['case_carriers_recovered'] and
      S['array_drivers']['case_wgs'] == _gc['case_carriers_wgs'], 'gene carrier totals differ from array.summary.json')
# the drivers typed directly on the arrays (05.rare_direct; not used by the gene readout, which is the imputed arrays
# only, user 2026-11-05; quoted in the notes when asked)
_rd = tsv('11.array', 'cohorts', 'full_mainland', '05.rare_direct', 'rare_direct.tsv')
S['array_typed_drivers'] = [{'variant': str(r.variant), 'array': r.array,
                             'group': r.group, 'n': int(r.n), 'wgs_carriers': int(r.wgs_carriers),
                             'both_carriers': int(r.both_carriers), 'array_only': int(r.array_only)}
                            for r in _rd.itertuples() if int(r.n) > 0 and (r.group == 'cases' or r.array == 'JSA24v1')]
_ir = [r['info_r2'] for r in _g if r['info_r2']]
S['array_lead_r2'] = [min(_ir), max(_ir)]
# where down-sampling acts (user 2026-11-05: "下取样的区域针对的是哪些位置"): the calling intervals, each signal
# variant and each tested gene variant +- 10 kb on the original 5-kb shard grid, with the signals each holds and how many
# signal variants fall in it; the HLA reads are those hla.typing extracted for HLA-HD (the classical MHC window, the
# chr6 ALT scaffolds, the HLA decoy contigs of hs38DH and the unmapped reads); nothing else is re-called
_hc = tsv('02.targets', 'hc_intervals.tsv')
_sv = tsv('02.targets', 'signal_variants.tsv')
_sv = _sv.drop_duplicates(['kind', 'unit', 'variant_id'])
S['ds_regions'] = []
for _i in _hc.itertuples():
    _in = _sv[(_sv.chrom == _i.chrom) & (_sv.pos >= _i.start) & (_sv.pos <= _i.end)]
    S['ds_regions'].append({'chrom': _i.chrom, 'start': int(_i.start), 'end': int(_i.end), 'kb': int(_i.length) // 1000,
                            'units': list(dict.fromkeys(_in.unit)),
                            'gwas_variants': int((_in.kind == 'gwas').sum()), 'gene_variants': int((_in.kind == 'gene').sum())})
check(sum(r['kb'] for r in S['ds_regions']) == S['design']['targets']['interval_kb'], 'interval kb differ from the design')
check(all(r['units'] for r in S['ds_regions']), 'a calling interval holds no signal variant')
_cl = [l.strip() for l in open(DF.ROOT.parent.parent / 'hla.typing' / 'results' / '00.manifest' / 'contig_list.txt')]
_mhc = [l for l in _cl if l.startswith('chr6:')]
check(len(_mhc) == 1, 'one MHC window expected in the HLA contig list')
_a, _b = _mhc[0].split(':')[1].split('-')
S['hla_reads'] = {'mhc': [int(_a), int(_b)], 'mhc_mb': [int(_a) / 1e6, int(_b) / 1e6], 'alt': sum(l.endswith('_alt') for l in _cl),
                  'decoy': sum(l.startswith('HLA') for l in _cl), 'unmapped': '*' in _cl, 'regions': len(_cl)}
S['array_counts'] = {c: {k: int((asig[asig.cohort == c].array_verdict == k).sum())
                         for k in asig.array_verdict.unique()} for c in CK}
S['array_contradicted'] = {c: int(asig[asig.cohort == c].array_verdict.str.contains('contradict').sum()) for c in CK}
S['array_consistent'] = {c: S['array_counts'][c].get('consistent on arrays', 0) for c in CK}
S['array_n_signals'] = {c: int((asig.cohort == c).sum()) for c in CK}

# ── control sources ─────────────────────────────────────────────────────────
cc = rtsv('controls.counts.tsv')
S['sources'] = {c: {r.source: r.controls_random for r in cc[(cc.cohort == c) & cc.subset.str.startswith('ctrlonly')].itertuples()}
                for c in CK}
csig = rtsv('controls.signals.tsv')
S['sources_signal'] = {}
for r in csig.itertuples():
    S['sources_signal'].setdefault(r.signal, {})[r.cohort] = {
        'p_Q': r.p_Q, 'I2': r.I2, 'or': r.or_pooled, 'weakest': r.weakest_loo_source, 'z': r.weakest_loo_z_obs,
        'z_exp': r.weakest_loo_z_expected, **{f'or_{s}': getattr(r, f'or_{s}') for s in SOURCES}}
S['sources_null'] = clean(rtsv('controls.null_lambda.tsv').to_dict('records'))
S['q_min']['sources'] = min(S['sources_signal'][g][c]['p_Q'] for g in qs for c in CK)

# the GEC threshold of the residue scan, as the discovery reading prints it, and to three digits
S['gec'] = {c: float(re.search(r'GEC ([\d.e-]+)', S['discovery']['DQB1:57:V'][c]['metric']).group(1)) for c in CK}
S['gec_exact'] = {c: S['ds_v57'][c]['p_orig'] / S['ds_v57'][c]['pg_orig'] for c in CK}
# the other thresholds of the main analyses, as their discovery readings print them: Bonferroni over the cohort's
# gene scan (per variant set); the residue positions are not on the slides
S['thresholds'] = {'genes': {}}
for g in S['discovery']:
    m = re.match(r'(\S+) \((\w+)\)$', g)
    if m:
        S['thresholds']['genes'].setdefault(m.group(2), {})
        for c in CK:
            S['thresholds']['genes'][m.group(2)][c] = float(re.search(r'Bonferroni ([\d.e-]+)',
                                                                      S['discovery'][g][c]['metric']).group(1))
# the number of ancestry components each cohort keeps (PopGMM: the cohorts are nested by these components)
_k = pd.read_csv(os.path.join(HERE, '..', '..', '..', 'PopGMM_output', 'keep_list_summary.tsv'), sep='\t')
S['cohort_k'] = {f'{r.variant}_mainland': int(r.n_components) for r in _k.itertuples() if f'{r.variant}_mainland' in CK}

# ── follow-up audits, every cohort ──────────────────────────────────────────
af = rtsv('followup.gwas_af.tsv')
S['lead_af'] = {}
for r in af.itertuples():
    S['lead_af'].setdefault(r.cohort, {}).setdefault(r.signal, {})[r.group] = {
        'n': None if pd.isna(r.n) else r.n, 'af': r.alt_freq}
gc = rtsv('followup.gene_carriers.tsv')
isgene = gc.unit.str.contains('(gene)', regex=False)
S['gene_carriers'], S['gene_sites'], S['panel_shares'] = {}, {}, {}
for r in gc[isgene].itertuples():
    S['gene_carriers'].setdefault(r.cohort, {}).setdefault(f'{r.gene}|{r.stratum}', {})[r.group] = {
        'n': r.n, 'carriers': r.carriers, 'share': r.carrier_share, 'hidden': r.hidden_as_missing,
        'tommo_expected': r.tommo_expected_carriers,
        'tommo_share': r.tommo_expected_carriers / r.n if r.n else None}
for r in gc[~isgene].itertuples():
    S['gene_sites'].setdefault(r.cohort, {}).setdefault(r.unit, {})[r.group] = {
        'n': r.n, 'carriers': r.carriers, 'share': r.carrier_share, 'hidden': r.hidden_as_missing,
        'tommo_expected': r.tommo_expected_carriers, 'tommo_af': r.tommo_af}
# the carrier share each reference panel predicts for a gene set: the sum over its tested variants of 2·AF (a variant
# the panel lacks counts as 0), beside the observed case and control shares
sites = gc[~isgene & (gc.group == 'case | all')]
for (c, g, st), d in sites.groupby(['cohort', 'gene', 'stratum']):
    k = f'{g}|{st}'
    x = {col: float((2 * pd.to_numeric(d[col], errors='coerce').fillna(0)).sum())
         for col in ('tommo_af', 'gnomad_eas_af', 'kg_eas_af', 'kg_jpt_af')}
    x.update(case=S['gene_carriers'][c][k]['case | all']['share'], control=S['gene_carriers'][c][k]['control | all']['share'],
             variants=int(len(d)))
    S['panel_shares'].setdefault(c, {})[k] = x
# the RPH3A driver (p.Arg544Cys): Nagahama against the other four sources, every cohort
S['rph3a_driver'] = {}
other = [f'control | {s}' for s in SOURCES if s != 'nagahama']
for c in CK:
    drv = S['gene_sites'][c]['chr12:112890858:C:T']
    d = {'nag_carriers': drv['control | nagahama']['carriers'], 'nag_n': drv['control | nagahama']['n'],
         'other_carriers': sum(drv[g]['carriers'] for g in other), 'other_n': sum(drv[g]['n'] for g in other),
         'ctrl_carriers': drv['control | all']['carriers'], 'ctrl_n': drv['control | all']['n'],
         'case_carriers': drv['case | all']['carriers'], 'case_n': drv['case | all']['n'],
         'tommo_expected': drv['control | all']['tommo_expected'], 'tommo_af': drv['control | all']['tommo_af'],
         'hidden_ctrl': drv['control | all']['hidden'],
         'hidden_by_source': {s: drv[f'control | {s}']['hidden'] for s in SOURCES}}
    for k in ('nag', 'other', 'ctrl', 'case'):
        d[f'{k}_share'] = d[f'{k}_carriers'] / d[f'{k}_n']
    d['tommo_share'] = d['tommo_expected'] / d['ctrl_n']
    S['rph3a_driver'][c] = d
# the largest number of possible carriers set missing at one site in one control source (all four genes)
srcs = [f'control | {s}' for s in SOURCES]
S['hidden_max_per_site_source'] = {c: int(max(d[g]['hidden'] for d in S['gene_sites'][c].values() for g in srcs if g in d))
                                   for c in CK}
# JAK2 p.Val617Phe
S['v617f'] = {}
for c in CK:
    x = S['gene_sites'][c]['chr9:5073770:G:T']
    S['v617f'][c] = {'case_carriers': x['case | all']['carriers'], 'case_hidden': x['case | all']['hidden'],
                     'case_reads': x['case | all']['carriers'] + x['case | all']['hidden'],
                     'ctrl_carriers': x['control | all']['carriers'], 'ctrl_hidden': x['control | all']['hidden'],
                     'tommo_expected': x['control | all']['tommo_expected'], 'case_n': x['case | all']['n'],
                     'ctrl_n': x['control | all']['n']}
gf = rtsv('followup.gene_filters.tsv')
S['gene_filters'] = {c: {'sites': int((gf.cohort == c).sum()), 'agree_all': bool((gf[gf.cohort == c].genotype_agreement == 1).all()),
                         'samples': int(gf[gf.cohort == c].samples.iloc[0]),
                         'missing_min': int(gf[gf.cohort == c].missing_published.min()),
                         'missing_max': int(gf[gf.cohort == c].missing_published.max())} for c in CK}
S['gene_filters_rows'] = int(len(gf))

# ── what each rule can detect, and the follow-up numbers the reviewers asked for (report.6 tables) ──────────────
from scipy import stats as SS  # noqa: E402


def allelic_z(a1, n1, a2, n2):
    """allelic test z: ALT copies a of N alleles in cases (1) and controls (2)"""
    p1, p2, p = a1 / n1, a2 / n2, (a1 + a2) / (n1 + n2)
    return (p1 - p2) / math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))


DEEP_GROUPS = ['case | NovaSeq 30x', 'case | DNBSeq-G400RS 30x']
# depth: the largest |z| drop Q1 could produce at a lead, if the deep cases' whole excess were an artifact (their
# ALT frequency set to the controls'); allelic test on followup.gwas_af.tsv (n called, ALT frequency)
S['depth_reach'] = {}
for c in CK:
    for g in gw:
        A = S['lead_af'][c][g]
        nc, ac = 2 * A['case | all']['n'], 2 * A['case | all']['n'] * A['case | all']['af']
        nk, ak = 2 * A['control | all']['n'], 2 * A['control | all']['n'] * A['control | all']['af']
        nd = sum(2 * A[x]['n'] for x in DEEP_GROUPS)
        ad = sum(2 * A[x]['n'] * A[x]['af'] for x in DEEP_GROUPS)
        z0 = allelic_z(ac, nc, ak, nk)
        z1 = allelic_z(ac - ad + nd * A['control | all']['af'], nc, ak, nk)
        # the DNBSeq-T7 cases' share of the cases' ALT excess over the controls: Q1 cannot test it (they sit at the
        # controls' depth already); Q2 and Q3 do (notes only, user 2026-11-04)
        nt7, at7 = 2 * A['case | DNBSeq-T7 30x']['n'], 2 * A['case | DNBSeq-T7 30x']['n'] * A['case | DNBSeq-T7 30x']['af']
        ex_all = ac - nc * A['control | all']['af']
        t7_share = (at7 - nt7 * A['control | all']['af']) / ex_all if ex_all else float('nan')
        S['depth_reach'].setdefault(c, {})[g] = {'z': abs(z0), 'z_null_deep': abs(z1), 'max_drop': abs(z0) - abs(z1),
                                                 'observed_drop': S['long'][g]['depth'][c]['value'],
                                                 't7_share': t7_share}
S['depth_reach_t7_share'] = {c: [min(x['t7_share'] for x in S['depth_reach'][c].values()),
                                 max(x['t7_share'] for x in S['depth_reach'][c].values())] for c in CK}
S['depth_reach_summary'] = {c: {'min': min(x['max_drop'] for x in S['depth_reach'][c].values()),
                                'max': max(x['max_drop'] for x in S['depth_reach'][c].values()),
                                'le1': int(sum(x['max_drop'] <= 1 for x in S['depth_reach'][c].values())),
                                'leads': len(S['depth_reach'][c])} for c in CK}
# each lead's observed drop as a share of its own reach (the readout the rule cannot give where the reach is <= 1),
# the leads where the rule (a drop > 1) can fire at all, and the leads whose replicate-median |z| does not drop
DEPTH_BAR = DF.CK.DEPTH_DZ                             # the rule's bar: a drop of more than 1 in |z|
for c in CK:
    xs = S['depth_reach'][c]
    share = {g: x['observed_drop'] / x['max_drop'] for g, x in xs.items() if x['max_drop'] > 0}
    for g, v in share.items():
        xs[g]['share'] = v
    top = max(share, key=share.get)
    S['depth_reach_summary'][c].update(
        max_share=share[top], max_share_lead=top,
        n_gt1=int(sum(x['max_drop'] > DEPTH_BAR for x in xs.values())),
        n_zero=int(sum(x['observed_drop'] <= 0 for x in xs.values())),
        n_negative=int(sum(x['observed_drop'] < 0 for x in xs.values())))
    check(S['depth_reach_summary'][c]['n_gt1'] + S['depth_reach_summary'][c]['le1'] == len(xs),
          f'{c}: depth reach counts do not add up')
S['t7_reads_lost'] = 1 - N['downsampled']['t7_fraction']

# the cohorts are nested: the cases and controls shared by all three, as shares of full mainland
stt = tsv('00.samples', 'sample_table.tsv', dtype=str)
inall = (stt.in_narrow_mainland == '1') & (stt.in_intermediate_mainland == '1') & (stt.in_full_mainland == '1')
S['cohort_overlap'] = {s: float((inall & (stt.status == s)).sum() / ((stt.in_full_mainland == '1') & (stt.status == s)).sum())
                       for s in ('case', 'control')}

# leave-one-out: the conditional model's spread and how many flags chance alone would give (evidence.loo*)
lo = loo.copy()
for col in ('z_pooled', 'sd_conditional', 'p_diff', 'testable', 'chance_of_flag'):
    check(col in lo.columns, f'evidence.loo.tsv lacks {col} (report.6)')
if 'sd_conditional' in lo.columns:
    lp = lo[(lo.check == 'lopo') & (lo.kind != 'gene')]
    S['loo_sd'] = {'t7': [float(lp[lp.left_out == 't7'].sd_conditional.min()), float(lp[lp.left_out == 't7'].sd_conditional.max())],
                   'small': [float(lp[lp.left_out != 't7'].sd_conditional.min()), float(lp[lp.left_out != 't7'].sd_conditional.max())]}
    # review round 9 (B-13): the same spread for the gene sets (CMC burden), which the range above leaves out
    _lgn = lo[(lo.check == 'lopo') & (lo.kind == 'gene') & lo.sd_conditional.notna()]
    S['loo_sd_genes'] = {'t7': [float(_lgn[_lgn.left_out == 't7'].sd_conditional.min()), float(_lgn[_lgn.left_out == 't7'].sd_conditional.max())],
                         'small': [float(_lgn[_lgn.left_out != 't7'].sd_conditional.min()), float(_lgn[_lgn.left_out != 't7'].sd_conditional.max())]}
cal = rtsv('evidence.loo_calibration.tsv')
S['loo_calibration'] = {f'{r.check}|{r.cohort}': {'testable': int(r.testable), 'observed': int(r.flags_observed),
                                                   'expected': float(r.flags_expected)} for r in cal.itertuples()}
S['loo_calibration_total'] = {'testable': int(cal.testable.sum()), 'observed': int(cal.flags_observed.sum()),
                              'expected': float(cal.flags_expected.sum())}
# per check (platform left out on the Q2 slide, source left out on the Q4 slide)
S['loo_calibration_by_check'] = {k: {'testable': int(x.testable.sum()), 'observed': int(x.flags_observed.sum()),
                                     'expected': float(x.flags_expected.sum())} for k, x in cal.groupby('check')}
# below the testability bar: a left-out group whose unchanged effect would still be significant (|z expected| at
# least 1.96 but under the bar of 2.8) and whose z falls under 1.96 (or turns): not read by the rule, disclosed once
zo, zb = DF.CK.LOO_Z_OBSERVED, DF.CK.LOO_Z_EXPECTED
bb = loo[(loo.z_expected.abs() >= zo) & (loo.z_expected.abs() < zb)
         & ((loo.z_observed.abs() < zo) | (loo.z_observed * loo.z_expected < 0))]
# the basis of each reading below 1.96: how likely so low a z would be if the effect were unchanged (z without the
# group ~ N(z_exp, sd_conditional), signed to the pooled effect); the same for the flag readings
from scipy.stats import norm as _N  # noqa: E402


def p_unchanged(r):
    sg = 1.0 if r.z_expected >= 0 else -1.0
    return float(_N.cdf((sg * r.z_observed - abs(r.z_expected)) / r.sd_conditional))


def power_cond(r):
    """P(|z| >= 1.96) for an unchanged effect, z ~ N(|z_exp|, sd_conditional): what the bar of 2.8 approximates
    with a spread of 1 (review round 7: the bar is conservative given the pooled estimate)"""
    m = abs(r.z_expected)
    return float(_N.sf((zo - m) / r.sd_conditional) + _N.cdf((-zo - m) / r.sd_conditional))


S['loo_below_bar'] = clean([{'signal': r.signal, 'kind': r.kind, 'check': r.check, 'cohort': r.cohort,
                             'left_out': r.left_out, 'z': r.z_observed, 'z_exp': r.z_expected, 'p_diff': r.p_diff,
                             'p_unchanged': p_unchanged(r), 'sd': r.sd_conditional, 'power_cond': power_cond(r)}
                            for r in bb.itertuples()])
fl = loo[(loo.z_expected.abs() >= zb) & ((loo.z_observed.abs() < zo) | (loo.z_observed * loo.z_expected < 0))]
S['loo_flag_basis'] = clean([{'signal': r.signal, 'check': r.check, 'cohort': r.cohort, 'left_out': r.left_out,
                              'z': r.z_observed, 'z_exp': r.z_expected, 'p_unchanged': p_unchanged(r), 'p_diff': r.p_diff}
                             for r in fl.itertuples()])
# gene sets left out (full mainland): the CMC burden z they are read with, and how many of their readings are testable
_gl = loo[(loo.kind == 'gene') & (loo.check == 'lopo') & (loo.cohort == 'full_mainland')]
S['gene_loo'] = {'readings': int(len(_gl)), 'not_testable': int((_gl.z_expected.abs() < zb).sum()),
                 'jak2_z': sorted(float(abs(v)) for v in _gl[_gl.signal.str.startswith('JAK2')].groupby('signal').z_pooled.first()),
                 't7_zexp_other': [float(_gl[(_gl.left_out == 't7') & ~_gl.signal.str.startswith('JAK2')].z_expected.abs().min()),
                                   float(_gl[(_gl.left_out == 't7') & ~_gl.signal.str.startswith('JAK2')].z_expected.abs().max())]}
# JAK2 M/H, the set called by SKAT-O, which the leave-one-out rule cannot read (no effect estimate): its SKAT-O P with
# each group left out and its case carriers by platform, described beside the rule (user 2026-11-06: a rationale for
# the rule as it is, not a new rule)
from scipy.stats import fisher_exact as _fe
S['jak2_skato'] = {}
for _c in CK:
    _j = {}
    for _sub in ['pooled', 'without_t7', 'without_novaseq', 'without_g400rs30', 'without_g400rs15', 'ctrlwithout_nagahama',
                 'ctrlwithout_bbj', 'ctrlwithout_htlv1', 'ctrlwithout_acc', 'ctrlwithout_lungca', 'only_t7', 'only_novaseq',
                 'only_g400rs30', 'only_g400rs15']:
        _d = pd.read_csv(os.path.join(R, '09.platform_assoc', _c, f'gene.{_sub}.moderate_high.tsv'), sep='\t')
        _d = _d[(_d.gene == 'JAK2') & (_d.arm == 'original')].set_index('method')
        _j[_sub] = {'p_skato': float(_d.loc['skato', 'pvalue']), 'p_cmc': float(_d.loc['cmc', 'pvalue']),
                    'case': int(_d.loc['skato', 'carriers_case']), 'control': int(_d.loc['skato', 'carriers_control']),
                    'rho': float(_d.loc['skato', 'rho'])}
    S['jak2_skato'][_c] = _j
_jf = S['jak2_skato']['full_mainland']
_cp = S['design']['case_platform']['full_mainland']
_t7c, _t7n = _jf['only_t7']['case'], int(_cp['DNBSeq-T7 30x'])
_oc = sum(_jf[f'only_{k}']['case'] for k in ('novaseq', 'g400rs30', 'g400rs15'))
_on = int(_cp['NovaSeq 30x'] + _cp['DNBSeq-G400RS 30x'] + _cp['DNBseq-G400RS 15x'])
S['jak2_platform_carriers'] = {'t7': _t7c, 't7_n': _t7n, 'other': _oc, 'other_n': _on,
                               't7_share': _t7c / _t7n, 'other_share': _oc / _on,
                               'fisher_p': float(_fe([[_t7c, _t7n - _t7c], [_oc, _on - _oc]])[1])}
check(_t7c + _oc == _jf['pooled']['case'], 'JAK2 case carriers by platform do not add up')
_gs = loo[(loo.kind == 'gene') & (loo.check == 'sources_loo') & (loo.cohort == 'full_mainland')]
S['gene_loo_src'] = {'readings': int(len(_gs)), 'not_testable': int((_gs.z_expected.abs() < zb).sum())}
S['loo_below_bar_n'] = {k: {c: int(((bb.check == k) & (bb.cohort == c)).sum()) for c in CK} for k in ('lopo', 'sources_loo')}
S['loo_below_bar_total'] = int(len(bb))
check(not bool(((bb.z_expected.abs() >= zb)).any()), 'loo_below_bar holds a testable group')

# readings: the grid, the defined readings, those that apply
# a check read once for the three cohorts (the gene sets' arrays: nested cohorts, counted in full mainland) is one
# merged cell where the others have three
_merged = int((chk.groupby(['signal', 'check']).size() == 1).sum())
S['grid'] = {'cells': S['n_signals'] * S['n_checks'] * len(CK) - 2 * _merged, 'defined': S['matrix']['readings'],
             'not_applicable': S['matrix']['cells'].get('not applicable', 0), 'merged': _merged}
S['grid']['apply'] = S['grid']['defined'] - S['grid']['not_applicable']
S['grid']['undefined'] = S['grid']['cells'] - S['grid']['defined']

# ── the evidence matrix as shown (user 2026-11-06: "现在的证据矩阵太让人误解了 ... 重新审视整个证据链条") ──────────
# The registry's seven checks are shown as seven readouts by question, with the same rules and flags:
#   - leaving out a case platform is split. One rule over four groups read as if the DNBSeq-T7 cases had been tested,
#     while that group (most of the cases) is testable for strong signals only and every other "consistent" came from
#     the three small platforms. So: without the DNBSeq-T7 cases (threat T3) and without one small platform (T2);
#   - Q and the random-effects estimate share one column: both compare the platforms' own estimates, common markers;
#   - a below-the-bar reading is marked in its cell.
# Readings per cell: consistent, flag, not testable, not checkable, not applicable, reported (gene carriers, one cell);
# 'below' marks a left-out group with 1.96 <= |z_exp| < 2.8 and |z| < 1.96. The union of the two left-out columns and
# the merged column reproduce the registry's readings (checked below).
M2_COLS = [('depth', 'Q1', ['T1', 'T4']), ('platform_est', 'Q2', ['T2', 'T3']), ('lopo_small', 'Q2', ['T2']),
           ('lopo_t7', 'Q2', ['T3']), ('array', 'Q3', ['T3', 'T4']), ('sources_Q', 'Q4', ['T5']),
           ('sources_loo', 'Q4', ['T5'])]
_lt = loo.assign(flagged=lambda d: d.testable.astype(bool) & ((d.z_observed.abs() < zo) | (d.z_observed * d.z_expected < 0)))
_bbk = {(r.signal, r.check, r.cohort, r.left_out) for r in bb.itertuples()}


def _left_out(sig, chk, c, groups):
    x = _lt[(_lt.signal == sig) & (_lt.check == chk) & (_lt.cohort == c) & _lt.left_out.isin(groups)]
    below = any((sig, chk, c, g) in _bbk for g in x.left_out)
    if not len(x):
        return {'r': 'not applicable', 'below': False}
    t = x[x.testable.astype(bool)]
    r = 'not testable' if not len(t) else ('flag' if t.flagged.any() else 'consistent')
    return {'r': r, 'below': below}


_SMALL = [k for k in SLUG if k != 't7']
_SRC = sorted(set(_lt[_lt.check == 'sources_loo'].left_out))
M2 = {'columns': [{'key': k, 'question': q, 'threats': t} for k, q, t in M2_COLS], 'rows': {}}
for g in [r['signal'] for r in S['signals']]:
    row = {}
    L_ = S['long'][g]

    def _rd(k, c, L_=L_):
        return L_.get(k, {}).get(c, {}).get('reading', 'not applicable')
    for k, _, _ in M2_COLS:
        if k == 'depth' or k == 'sources_Q':
            row[k] = {c: {'r': _rd(k, c), 'below': False} for c in CK}
        elif k == 'platform_est':
            cell = {}
            for c in CK:
                q_, re_ = _rd('platform_Q', c), _rd('platform_RE', c)
                check((q_ == 'not applicable') == (re_ == 'not applicable'), f'matrix2: {g} {c}: Q and RE apply apart')
                cell[c] = {'r': 'flag' if 'flag' in (q_, re_) else q_, 'below': False}
            row[k] = cell
        elif k == 'lopo_small':
            row[k] = {c: _left_out(g, 'lopo', c, _SMALL) for c in CK}
        elif k == 'lopo_t7':
            row[k] = {c: _left_out(g, 'lopo', c, ['t7']) for c in CK}
        elif k == 'sources_loo':
            row[k] = {c: _left_out(g, 'sources_loo', c, _SRC) for c in CK}
        elif k == 'array':
            a = L_.get('array', {})
            if len(a) == 1:                     # the gene sets: one cell for the three nested array subsets
                (c1, v), = a.items()
                m_ = re.search(r'(\d+) of (\d+)', str(v['metric']))
                row[k] = {'merged': {'r': v['reading'], 'below': False,
                                     'count': [int(m_.group(1)), int(m_.group(2))] if m_ else None}}
            else:
                row[k] = {c: {'r': {'consistent on arrays': 'consistent', 'contradicted': 'flag'}.get(
                    a[c]['reading'], a[c]['reading']), 'below': False} for c in CK}
    # the split and merged columns reproduce the registry's readings
    for c in CK:
        s_, t_ = row['lopo_small'][c]['r'], row['lopo_t7'][c]['r']
        u_ = 'flag' if 'flag' in (s_, t_) else ('not testable' if s_ == t_ == 'not testable' else 'consistent')
        check(u_ == _rd('lopo', c), f'matrix2: {g} {c}: the split lopo ({s_}, {t_}) does not give {_rd("lopo", c)}')
        check(row['sources_loo'][c]['r'] == _rd('sources_loo', c), f'matrix2: {g} {c}: sources_loo differs')
    M2['rows'][g] = row
# counts: cells by reading (a merged cell counts once), the below-the-bar marks, per column the signals the rule reads
# (consistent or flag) in each cohort and the flag cells
_cells = [v for row in M2['rows'].values() for k, cell in row.items() for v in cell.values()]
M2['cells'] = {}
for v in _cells:
    M2['cells'][v['r']] = M2['cells'].get(v['r'], 0) + 1
M2['n_cells'] = len(_cells)
M2['below'] = sum(v['below'] for v in _cells)
M2['read'] = {k: {c: sum(M2['rows'][g][k].get(c, {}).get('r') in ('consistent', 'flag') for g in M2['rows']) for c in CK}
              for k, _, _ in M2_COLS}
M2['flags'] = {k: sum(v['r'] == 'flag' for g in M2['rows'] for v in M2['rows'][g][k].values()) for k, _, _ in M2_COLS}
check(M2['n_cells'] == S['grid']['cells'], f'matrix2: {M2["n_cells"]} cells against the grid\'s {S["grid"]["cells"]}')
check(M2['below'] == len({(r.signal, r.check, r.cohort, 'lopo_t7' if (r.check == 'lopo' and r.left_out == 't7') else r.check)
                          for r in bb.itertuples()}), 'matrix2: the below-the-bar marks do not match loo_below_bar')
check(sum(M2['flags'].values()) == S['flag_cells'], 'matrix2: the flag cells differ from the registry\'s')
S['matrix2'] = M2

# arrays in the matrix's units: signals with an array reading, and how many it supports, per cohort
arr = lg[(lg.check == 'array') & (lg.reading != 'not applicable')]
S['array_matrix'] = {c: {'covered': int(((arr.cohort == c) & (arr.kind != 'gene')).sum()),   # association tests
                         'consistent': int(((arr.cohort == c) & (arr.reading == 'consistent') & (arr.kind != 'gene')).sum()),
                         'flag': int(((arr.cohort == c) & (arr.reading == 'flag')).sum())} for c in CK}
# ... and how many of them the arrays could contradict at all (scripts/array_evidence.py): an association row (typed,
# imputed, LD proxy, HIBAG) can be contradicted; no evidence (an indel without a proxy) cannot; the gene sets are
# reported, not ruled on (user 2026-11-05: their driver carriers are counted on the imputed arrays, full mainland)
ASSOC_EV = ('typed', 'imputed', 'imputed (HIBAG)', 'LD proxy')
kind_of = {r['signal']: r['kind'] for r in S['signals']}
for c in CK:
    n_assoc = 0
    for sig in arr[arr.cohort == c].signal:
        e = S['array_signal'][c].get(sig.split(' (')[0], {}).get('evidence')
        n_assoc += kind_of[sig] != 'gene' and e in ASSOC_EV
    S['array_matrix'][c]['contradictable'] = int(n_assoc)
    S['array_matrix'][c]['reported_genes'] = int(((arr.cohort == c) & (arr.reading == 'reported')).sum())
# the arrays' verdicts in the matrix's units (gene sets apart), every cohort
amap = {r.signal: r.reading for r in arr.itertuples()}
for c in CK:
    a_ = arr[arr.cohort == c]
    units = {}
    for sig in a_.signal:
        units[sig] = S['array_signal'][c][sig.split(' (')[0]]['verdict']
    S['array_matrix'][c]['verdicts'] = {v: int(sum(1 for x in units.values() if x == v)) for v in sorted(set(units.values()))}

# the signals each question can rule on: at least one of its checks reads "consistent" or "flag" in some cohort
# (JAK2's gene sets are not testable when groups are left out; position 57 is reported, without a rule)
QCHECKS = {'Q2': ['platform_Q', 'platform_RE', 'lopo'], 'Q4': ['sources_Q', 'sources_loo']}
S['q_rulable'] = {}
for q, ks in QCHECKS.items():
    sub = chk[chk.check.isin(ks)]
    rul = sorted(s for s, x in sub.groupby('signal') if x.reading.isin(['consistent', 'flag']).any())
    flg = sorted(r.signal for r in ev.itertuples() if set(re.split(r'[;,\s]+', str(r.flags))) & set(ks))
    S['q_rulable'][q] = {'signals': len(rul), 'flagged': len(flg), 'not_flagged': len(set(rul) - set(flg)),
                         'unreadable': sorted(set(ev.signal) - set(rul)), 'flagged_signals': flg}
    check(set(flg) <= set(rul), f'{q}: a flagged signal the question cannot rule on')

# heterogeneity: INPP5D, the rare lead outside the Q rule (pooled ALT frequency, empirical Q P, case-only fit)
ps = tsv('10.platform', '_summary', 'platform.signals.tsv')
S['inpp5d'] = {}
for c in CK:
    A = S['lead_af'][c]['INPP5D']
    na, nk = A['case | all']['n'], A['control | all']['n']
    co = tsv('10.platform', c, 'platform.case_only.tsv')
    r = co[co.marker == S['lead_marker']['INPP5D']].iloc[0]
    q = ps[(ps.cohort == c) & (ps.signal == 'INPP5D')].iloc[0]
    S['inpp5d'][c] = {'maf': (na * A['case | all']['af'] + nk * A['control | all']['af']) / (na + nk),
                      'p_Q_emp': float(q.p_Q_empirical), 'p_case_only': None if pd.isna(r.p_lrt) else float(r.p_lrt),
                      'converged': str(r.converged) == 'True' and not pd.isna(r.p_lrt)}
# why heterogeneity is read at MAF >= 5% only: a rare lead's case ALT copies in the three small platform strata, and
# the one-step stratum OR they give (full mainland; platform.strata.tsv, original arm)
_st = tsv('10.platform', 'full_mainland', 'platform.strata.tsv')
_st = _st[_st.arm == 'original']
_lm = {v: k for k, v in S['lead_marker'].items()}
_st = _st[_st.variant_id.isin(_lm)].assign(lead=lambda d: d.variant_id.map(_lm))
_st['alt_case'] = (2 * _st.n_case * _st.af_case).round().astype(int)
_ctrl_af = {g: S["lead_af"]["full_mainland"][g]["control | all"]["af"] for g in _lm.values() if g in S["lead_af"]["full_mainland"]}
_rare = [g for g, f in _ctrl_af.items() if min(f, 1 - f) < DF.CK.Q_MIN_MAF]
_small = _st[_st.lead.isin(_rare) & _st.subset.isin(['only_novaseq', 'only_g400rs30', 'only_g400rs15'])]
_ip = _st[(_st.lead == 'INPP5D') & (_st.subset == 'only_g400rs30')].iloc[0]
S['rare_small_strata'] = {'leads': len(_rare), 'alt_min': int(_small.alt_case.min()), 'alt_max': int(_small.alt_case.max()),
                          'inpp5d_g400rs30_alt': int(_ip.alt_case), 'inpp5d_g400rs30_or': float(_ip['or'])}
# the two layers side by side: a common lead (FTO) and a rare one (INPP5D) in the three small strata: copies and OR
S['strata_example'] = {}
for _g in ('FTO', 'INPP5D'):
    _x = _st[(_st.lead == _g) & _st.subset.isin(['only_novaseq', 'only_g400rs30', 'only_g400rs15'])]
    S['strata_example'][_g] = {'alt': [int(_x.alt_case.min()), int(_x.alt_case.max())],
                               'or': [float(_x['or'].min()), float(_x['or'].max())],
                               'maf': float(S['lead_af']['full_mainland'][_g]['control | all']['af'])}
# Q4's two layers: the same two leads in the control sources (every case shared; full mainland, original arm):
# control ALT copies and the one-step OR of each source alone; and the rare leads' copies in the three smallest sources
_cs = tsv('12.controls', 'controls.strata.tsv')
_cs = _cs[(_cs.arm == 'original') & (_cs.cohort == 'full_mainland') & _cs.subset.str.startswith('ctrlonly_')]
_cs = _cs[_cs.variant_id.isin(_lm)].assign(lead=lambda d: d.variant_id.map(_lm))
_cs['alt_ctrl'] = (2 * _cs.n_ctrl * _cs.af_ctrl).round().astype(int)
S['source_strata_example'] = {}
for _g in ('FTO', 'INPP5D'):
    _x = _cs[_cs.lead == _g]
    check(len(_x) == 5, f'{_g}: five control sources expected')
    S['source_strata_example'][_g] = {'alt': [int(_x.alt_ctrl.min()), int(_x.alt_ctrl.max())],
                                      'or': [float(_x['or'].min()), float(_x['or'].max())],
                                      'maf': float(S['lead_af']['full_mainland'][_g]['control | all']['af'])}
_small3 = sorted(S['sources']['full_mainland'], key=lambda k: S['sources']['full_mainland'][k])[:3]
_xs = _cs[_cs.lead.isin(_rare) & _cs.subset.isin([f'ctrlonly_{k}' for k in _small3])]
S['source_rare_small'] = {'leads': len(_rare), 'sources': _small3,
                          'n': [int(min(S['sources']['full_mainland'][k] for k in _small3)),
                                int(max(S['sources']['full_mainland'][k] for k in _small3))],
                          'alt_min': int(_xs.alt_ctrl.min()), 'alt_max': int(_xs.alt_ctrl.max())}
S['case_only_nonconverged'] = {c: int((tsv('10.platform', c, 'platform.case_only.tsv').converged.astype(str) != 'True').sum())
                               for c in CK}

# reference panel: the control deficit against ToMMo's prediction (Poisson, lower tail), every gene set
S['ctrl_deficit'] = {}
for c in CK:
    for k, d in S['gene_carriers'][c].items():
        x = d['control | all']
        S['ctrl_deficit'].setdefault(k, {})[c] = {'obs': int(x['carriers']), 'exp': float(x['tommo_expected']),
                                                  'p': float(SS.poisson.cdf(x['carriers'], x['tommo_expected'])),
                                                  'fold': float(x['tommo_expected'] / x['carriers']) if x['carriers'] else None}

# RPH3A driver: Nagahama against the other sources (Fisher) and against ToMMo's prediction (Poisson)
for c in CK:
    d = S['rph3a_driver'][c]
    d['fisher_p'] = float(SS.fisher_exact([[d['nag_carriers'], d['nag_n'] - d['nag_carriers']],
                                           [d['other_carriers'], d['other_n'] - d['other_carriers']]])[1])
    d['nag_expected'] = d['tommo_share'] * d['nag_n']
    d['nag_poisson_p'] = float(SS.poisson.cdf(d['nag_carriers'], d['nag_expected']))

# ADARB2: exact tests on allele counts (DNBSeq-T7 cases against the other case groups and the controls), the
# pre-QC records of the DNBSeq-T7 cases, and the records within 50 bp of the lead
OTHER_CASE = ['case | NovaSeq 30x', 'case | DNBSeq-G400RS 30x', 'case | DNBseq-G400RS 15x']
S['adarb2'] = {}
lpq = rtsv('followup.lead_preqc.tsv')
lw = rtsv('followup.lead_window.tsv')
for c in CK:
    A = S['lead_af'][c]['ADARB2']
    al = lambda gs: (round(sum(2 * A[x]['n'] * A[x]['af'] for x in gs)), round(sum(2 * A[x]['n'] for x in gs)))  # noqa: E731
    t7, oc, ct = al(['case | DNBSeq-T7 30x']), al(OTHER_CASE), al(['control | all'])
    pq = lpq[(lpq.cohort == c) & (lpq.signal == 'ADARB2') & (lpq.group == 'case | DNBSeq-T7 30x')].iloc[0]
    w = lw[(lw.cohort == c) & (lw.signal == 'ADARB2')]
    lead = w[w.is_lead.astype(str) == 'True'].set_index('group').carriers
    d28 = w[w.variant_id.str.startswith('chr10:1747171:')].set_index('group').carriers
    snv = w[(w.offset_bp.abs() <= 50) & (w.variant_id.str.split(':').str[2].str.len() == 1)
            & (w.variant_id.str.split(':').str[3].str.len() == 1)]
    S['adarb2'][c] = {'t7_alt': t7[0], 't7_an': t7[1], 'other_alt': oc[0], 'other_an': oc[1], 'ctrl_alt': ct[0], 'ctrl_an': ct[1],
                      'p_t7_other': float(SS.fisher_exact([[t7[0], t7[1] - t7[0]], [oc[0], oc[1] - oc[0]]])[1]),
                      'p_t7_ctrl': float(SS.fisher_exact([[t7[0], t7[1] - t7[0]], [ct[0], ct[1] - ct[0]]])[1]),
                      't7_pre_carriers': int(pq.pre_carriers), 't7_alt_read_samples': int(pq.alt_read_samples),
                      't7_n': int(pq.n),
                      # review round 9 (B-05): the genotype filters at the lead (the record itself passes the site QC)
                      'qc_missing': {g: {'n': int(r.set_missing_by_filters), 'alt': int(r.set_missing_with_alt),
                                         'pre_carriers': int(r.pre_carriers)}
                                     for g, r in lpq[(lpq.cohort == c) & (lpq.signal == 'ADARB2')].set_index('group').iterrows()
                                     if g in ('case | all', 'control | all', 'case | DNBSeq-T7 30x')},
                      'del28_same_counts': bool(lead.reindex(d28.index).eq(d28).all() and len(d28) == len(lead)),
                      'del28_ctrl': int(d28['control | all']), 'del28_t7': int(d28['case | DNBSeq-T7 30x']),
                      'del28_bp': len(w[w.variant_id.str.startswith('chr10:1747171:')].variant_id.iloc[0].split(':')[2]) - 1,
                      'snv_t7': {int(r.offset_bp): int(r.carriers) for r in snv[snv.group == 'case | DNBSeq-T7 30x'].itertuples()
                                 if r.carriers > 0},
                      # review round 8 (B-11): the same SNVs' carrier shares in the T7 cases and in the controls
                      'snv_t7_pct': {int(r.offset_bp): float(r.carriers / r.n)
                                     for r in snv[snv.group == 'case | DNBSeq-T7 30x'].itertuples() if r.carriers > 0},
                      'snv_ctrl_pct': {int(r.offset_bp): float(r.carriers / r.n)
                                       for r in snv[snv.group == 'control | all'].itertuples() if r.carriers > 0}}
# indel inflation of the DNBSeq-T7 stratum on the null markers (mean over the five MAF bins)
S['t7_null_lambda'] = {}
for c in CK:
    nc_ = tsv('10.platform', c, 'platform.null_calibration.tsv')
    t = nc_[nc_.subset == 'only_t7']
    S['t7_null_lambda'][c] = {'indel': float(t[t.is_indel.astype(str) == 'True']['lambda'].mean()),
                              'snv': float(t[t.is_indel.astype(str) == 'False']['lambda'].mean())}
gv = exp[(exp.kind == 'gwas') & (exp.cohort == 'full_mainland')].variant_id
S['gwas_indels'] = int(sum(len(v.split(':')[2]) != 1 or len(v.split(':')[3]) != 1 for v in gv))

# missing genotypes at the leads: cases against controls (published genotypes)
mis = lpq[lpq.group.isin(['case | all', 'control | all']) & (lpq.signal.isin(gw))].copy()
mis['share'] = mis.published_missing / mis.n
S['lead_missing'] = {c: {g.split(' | ')[0]: float(mis[(mis.cohort == c) & (mis.group == g)].share.median())
                         for g in ('case | all', 'control | all')} for c in CK}

# LINC01146 NovaSeq: the standard error of its ALT frequency (n = 50)
for c in CK:
    A = S['lead_af'][c]['LINC01146']['case | NovaSeq 30x']
    A['se'] = math.sqrt(A['af'] * (1 - A['af']) / (2 * A['n']))
S['platform_diff_tests'] = {c: int(sum(ps[(ps.cohort == c) & (ps.role == 'Lead')][f'p_diff_{k}'].notna().sum() for k in SLUG))
                             for c in CK}

# control sources: a source that differs from the rest at a lead (P of the difference between that source alone and
# the controls without it < 0.001, a descriptive cut, not a rule), and in which direction: it weakens the signal when
# its own estimate lies between 0 and the estimate without it (or on the other side of 0), it strengthens the signal
# when its estimate lies further from 0 than the estimate without it
cl = tsv('12.controls', 'controls.loo.tsv')
m2s = {v: k for k, v in S['lead_marker'].items()}
DIFF_P = 0.001
dif = cl[cl.variant_id.isin(m2s) & (cl.p_diff < DIFF_P)]
S['source_differs'] = [{'lead': m2s[r.variant_id], 'source': r.left_out, 'cohort': r.cohort, 'p_diff': float(r.p_diff),
                        'beta_source': float(r.beta_source_only), 'beta_rest': float(r.beta),
                        'direction': 'strengthens' if (r.beta_source_only * r.beta > 0
                                                       and abs(r.beta_source_only) > abs(r.beta)) else 'weakens'}
                       for r in dif.itertuples()]
S['source_differs_n'] = {c: {d: len({x['lead'] for x in S['source_differs'] if x['cohort'] == c and x['direction'] == d})
                             for d in ('weakens', 'strengthens')} for c in CK}
S['source_differs_pairs'] = sorted({(x['lead'], x['source'], x['direction']) for x in S['source_differs']})
S['source_differs_cut'] = DIFF_P
# case platforms: one platform against the rest at a lead (P_diff of the leave-one-out table), descriptive, as for the
# sources (review round 7: two rare leads differ in a small platform; the rule does not read P_diff)
_lp = loo[(loo.check == 'lopo') & (loo.kind == 'gwas') & loo.p_diff.notna()]
S['platform_pdiff'] = {f'{r.signal}|{r.left_out}|{r.cohort}': float(r.p_diff) for r in _lp.itertuples()}
S['platform_differs'] = sorted({(r.signal, r.left_out) for r in _lp[_lp.p_diff < DIFF_P].itertuples()})
# review round 7, U4 (2026-11-06): at a rare lead (MAF < 5%, outside Q) a small platform's one-step estimate overstates
# P_diff; the descriptive counterpart counts that platform's ALT alleles against the other cases' (Fisher's exact test on
# the called genotypes of followup.gwas_af.tsv), for every rare lead, case platform and cohort. Not a rule.
S['rare_leads'] = [r['signal'] for r in S['signals'] if r['kind'] == 'gwas' and r['signal'] not in S['q_markers']]


def _alt_alleles(c, g, grp):
    x = S['lead_af'][c][g][grp]
    n2 = 2 * x['n']
    a = x['af'] * n2
    check(abs(a - round(a)) < 1e-6, f'rare_platform_exact: {g} {grp} {c}: the ALT count {a} is not an integer')
    return int(round(a)), int(round(n2))


S['rare_platform_exact'] = {}
for _g in S['rare_leads']:
    for c in CK:
        _A, _N2 = _alt_alleles(c, _g, 'case | all')
        _parts = {k: _alt_alleles(c, _g, f'case | {p}') for k, p in SLUG.items()}
        check(sum(v[0] for v in _parts.values()) == _A and sum(v[1] for v in _parts.values()) == _N2,
              f'rare_platform_exact: {_g} {c}: the case platforms do not add up to all cases')
        for _k, (_a, _n2) in _parts.items():
            S['rare_platform_exact'][f'{_g}|{_k}|{c}'] = {
                'alt': _a, 'alleles': _n2, 'alt_rest': _A - _a, 'alleles_rest': _N2 - _n2,
                'p': float(SS.fisher_exact([[_a, _n2 - _a], [_A - _a, (_N2 - _n2) - (_A - _a)]])[1])}
# the pairs (lead, platform) below 0.05 in some cohort, and how many pairs were tested per cohort
S['rare_platform_exact_lt05'] = sorted({tuple(k.split('|')[:2]) for k, v in S['rare_platform_exact'].items()
                                        if v['p'] < 0.05})
S['rare_platform_exact_pairs'] = len(S['rare_leads']) * len(SLUG)
check(all(f'{g}|{k}|{c}' in S['rare_platform_exact'] for g, k in S['platform_differs'] for c in CK),
      'rare_platform_exact: a platform that differs at a lead is not a rare lead')
S['source_differs_p'] = {f'{m2s[r.variant_id]}|{r.left_out}|{r.cohort}': float(r.p_diff) for r in
                         cl[cl.variant_id.isin(m2s)].itertuples()}
check(all(len({x['direction'] for x in S['source_differs'] if (x['lead'], x['source']) == (g, s)}) == 1
          for g, s, _ in S['source_differs_pairs']), 'a source changes direction between cohorts')

# Q4's own logic (user 2026-11-04: "Q4也有其特殊性，不要一味一样"): every case stays in every subset, so leaving out a
# source costs little power; per source what leaving it out costs (median z_exp / z_pooled over the GWAS leads and HLA
# markers, range over cohorts) and finds (testable readings, flags, the chance of a flag, three cohorts summed);
# the same ratio without the DNBSeq-T7 cases for comparison
_sl = loo[loo.check == 'sources_loo'].copy()
_sl['ratio'] = _sl.z_expected.abs() / _sl.z_pooled.abs()
_sl['is_flag'] = (_sl.z_expected.abs() >= zb) & ((_sl.z_observed.abs() < zo) | (_sl.z_observed * _sl.z_expected < 0))
S['src_loo_power'] = {}
for _k in S['sources']['full_mainland']:
    _a = _sl[_sl.left_out == _k]
    _r = _a[_a.kind != 'gene'].groupby('cohort').ratio.median()
    _t = _a[_a.z_expected.abs() >= zb]
    S['src_loo_power'][_k] = {'ratio': [float(_r.min()), float(_r.max())], 'testable': int(len(_t)),
                              'expected': float(_t.chance_of_flag.sum()), 'observed': int(_a.is_flag.sum())}
check(abs(sum(v['expected'] for v in S['src_loo_power'].values()) - S['loo_calibration_by_check']['sources_loo']['expected']) < 1e-6,
      'per-source chances do not sum to the calibration total')
# every reading of leaving a source out, by signal type (three cohorts): the slide's point is that none changes
_bbk = (_sl.z_expected.abs() < zb) & (_sl.z_expected.abs() >= zo) & (_sl.z_observed.abs() < zo)
S['src_loo_kinds'] = {}
for _kd, _x in _sl.groupby('kind'):
    _nt = _x.z_expected.abs() < zb
    S['src_loo_kinds'][_kd] = {'signals': int(_x.signal.nunique()), 'readings': int(len(_x)),
                               'consistent': int((~_nt & ~_x.is_flag).sum()), 'not_testable': int(_nt.sum()),
                               'below_bar': int((_bbk & (_sl.kind == _kd)).sum()), 'flag': int(_x.is_flag.sum())}
check(sum(v['below_bar'] for v in S['src_loo_kinds'].values()) ==
      sum(1 for r in S['loo_below_bar'] if r['check'] == 'sources_loo'), 'below-bar count differs from loo_below_bar')
_t7 = loo[(loo.check == 'lopo') & (loo.left_out == 't7') & (loo.kind != 'gene')]
_t7r = (_t7.z_expected.abs() / _t7.z_pooled.abs()).groupby(_t7.cohort).median()
S['t7_loo_ratio'] = [float(_t7r.min()), float(_t7r.max())]
# leaving the DNBSeq-T7 cases out: the signals for which that reading is testable (review round 7: the Q2 answer and
# the matrix notes are qualified by it)
_t7a = loo[(loo.check == 'lopo') & (loo.left_out == 't7') & (loo.kind != 'hla_position')]
S['t7_loo_testable'] = {c: {'signals': int(_t7a[(_t7a.cohort == c) & (_t7a.z_expected.abs() >= zb)].signal.nunique()),
                            'gwas': int(_t7a[(_t7a.cohort == c) & (_t7a.kind == 'gwas')
                                             & (_t7a.z_expected.abs() >= zb)].signal.nunique()),
                            'of': int(_t7a[_t7a.cohort == c].signal.nunique())} for c in CK}
check(all(S['t7_loo_testable'][c]['of'] == S['n_signals'] for c in CK), 't7_loo_testable: a signal lacks its T7 reading')
check(S['matrix2']['read']['lopo_t7'] == {c: S['t7_loo_testable'][c]['signals'] for c in CK},
      'matrix2: the signals read without the DNBSeq-T7 cases differ from t7_loo_testable')
# gene sets left out of a source (full mainland): nt because the pooled burden is weak, or because of the left-out source
_gf = _sl[(_sl.kind == 'gene') & (_sl.cohort == 'full_mainland')]
S['gene_loo_src'].update({
    'weak_pooled': int((_gf.z_pooled.abs() < zb).sum()),
    'weak_sets': sorted(set(_gf[_gf.z_pooled.abs() < zb].signal)),
    'by_loss': int(((_gf.z_pooled.abs() >= zb) & (_gf.z_expected.abs() < zb)).sum()),
    'osbp_mh_nagahama': {c: float(abs(_sl[(_sl.signal == 'OSBP (moderate_high)') & (_sl.cohort == c) &
                                          (_sl.left_out == 'nagahama')].z_expected.iloc[0])) for c in CK}})
check(S['gene_loo_src']['weak_pooled'] + S['gene_loo_src']['by_loss'] == S['gene_loo_src']['not_testable'],
      'gene nt readings are neither weak pooled nor lost by leaving out')
# each source that differs from the rest at a lead: both estimates, the cohorts where it differs, the reading without it
S['source_pairs'] = []
for _g, _k, _d in S['source_differs_pairs']:
    _v = S['lead_marker'][_g]
    _row = {'lead': _g, 'source': _k, 'direction': _d,
            'cohorts': [c for c in CK if S['source_differs_p'][f'{_g}|{_k}|{c}'] < DIFF_P], 'by_cohort': {}}
    for c in CK:
        _r = cl[(cl.variant_id == _v) & (cl.left_out == _k) & (cl.cohort == c)].iloc[0]
        _e = _sl[(_sl.signal == _g) & (_sl.cohort == c) & (_sl.left_out == _k)].iloc[0]
        _sg = 1.0 if _e.z_expected >= 0 else -1.0
        _row['by_cohort'][c] = {'p_diff': float(_r.p_diff), 'or_source': math.exp(_r.beta_source_only),
                                'or_rest': math.exp(_r.beta), 'z': float(_sg * _e.z_observed), 'z_exp': float(abs(_e.z_expected))}
    S['source_pairs'].append(_row)

# imputed leads not confirmed person by person (case r²) and their frequencies, every cohort
S['imputed_leads'] = {}
for c in CK:
    il = tsv('11.array', 'cohorts', c, '06.imputed', 'leads.tsv')
    for r in il.itertuples():
        S['imputed_leads'].setdefault(r.signal, {})[c] = {
            'r2_cases': r.r2_cases, 'case_wgs': r.freq_cases_wgs, 'case_imp': r.freq_cases_imputed,
            'ctrl_wgs': r.freq_controls_wgs, 'ctrl_imp': r.freq_controls_imputed}
S['array_dz4'] = {c: int(N['array'][c]['abs_dz_gt4']) for c in CK}

# JAK2 p.Val617Phe: carriers in the DNBSeq-T7 cases against the other case groups (Fisher)
for c in CK:
    x = S['gene_sites'][c]['chr9:5073770:G:T']
    t7c, t7n = x['case | DNBSeq-T7 30x']['carriers'], x['case | DNBSeq-T7 30x']['n']
    oc_ = sum(x[g]['carriers'] for g in OTHER_CASE)
    on_ = sum(x[g]['n'] for g in OTHER_CASE)
    S['v617f'][c].update(t7_carriers=int(t7c), t7_n=int(t7n), other_carriers=int(oc_), other_n=int(on_),
                         p_t7_other=float(SS.fisher_exact([[t7c, t7n - t7c], [oc_, on_ - oc_]])[1]))

# the platform strata's shared-control correlation: analytic against the null-SNP correlation (DNBSeq-T7 pairs)
rh = pd.concat([tsv('10.platform', c, 'platform.rho.tsv') for c in CK])
rt = rh[(rh.stratum_1 == 'only_t7') | (rh.stratum_2 == 'only_t7')]
S['rho_t7'] = {'analytic': [float(rt.rho_analytic.min()), float(rt.rho_analytic.max())],
               'null': [float(rt.r_null_snps.min()), float(rt.r_null_snps.max())]}

# review round 7, U5 (2026-11-06): the calibration of the two control-source checks over the null SNPs. The pipeline
# reads the platforms' Q against the null SNPs (platform_stats.py) but the sources' Q by its analytic chi-square P, and it
# checks the leave-one-out model on the null SNPs for neither; both are computed here from the subsets' SAIGE rows
# (09.platform_assoc) with the pipeline's own pieces (platform_stats.rho_shared, maf_bin and the GLS Q of
# control_report.het_rows, validated below against controls.heterogeneity.tsv). Descriptive (notes and one caption); the
# rules and their readings are unchanged.
#   sources_null_q   per cohort: null SNPs at MAF >= 5% (the Q rule's bar) and the share with analytic source-Q P < 0.05
#                    and < the threshold
#   sources_q_emp    per cohort x common marker: the analytic source-Q P and the P read against the null SNPs of the
#                    same MAF bin and type, (1 + #(Q_null >= Q)) / (1 + n), as platform_stats.py does for the platforms
#   loo_null_sd      SD over the null SNPs of (z - z_exp) / sd_conditional, the leave-one-out model's standardized
#                    deviation (1 if the model holds), across the left-out platforms, the sources, and Nagahama alone
#   loo_null_undefined  the largest share of null SNPs whose SE without a source is not above the pooled SE (so
#                    sd_conditional is undefined), with its source and cohort
import numpy as np  # noqa: E402
import platform_stats as PS  # noqa: E402  (the pipeline's formulas; read only, no bytecode written)
import report_lib as RL  # noqa: E402
_nul = pd.read_csv(os.path.join(R, '02.targets', 'null_snps.tsv'), sep='\t', dtype={'variant_id': str}).set_index('variant_id')
_hq = tsv('12.controls', 'controls.heterogeneity.tsv', dtype={'variant_id': str})
_COMMON = {f'{lo}-{hi}' for lo, hi in zip(PS.MAF_BINS[:-1], PS.MAF_BINS[1:]) if lo >= DF.CK.Q_MIN_MAF}
S['sources_null_q'], S['sources_q_emp'] = {}, {}
_sd, _und = {'platforms': [], 'sources': [], 'nagahama': []}, []
for c in CK:
    _cnt = pd.read_csv(RL.subset_counts_file(R, c, RL.SOURCE_ONLY[0]), sep='\t').set_index('subset')
    _pc = pd.read_csv(RL.subset_counts_file(R, c, 'pooled'), sep='\t').set_index('subset')
    _n = int(_pc.loc['pooled', 'cases_random'])
    _m = {s: int(_cnt.loc[f'ctrlonly_{s}', 'controls_random']) for s in RL.SOURCES}
    _sa = RL.read_saige_subsets(R, c, ['pooled'] + RL.SOURCE_ONLY + RL.SOURCE_WITHOUT + RL.PLATFORM_WITHOUT)
    check(len(_sa) == 1 + 2 * len(RL.SOURCES) + len(RL.PLATFORM_WITHOUT), f'sources_null_q: {c}: a subset is missing')
    _Ri = np.linalg.inv(np.array([[1.0 if i == j else PS.rho_shared(_m[a], _m[b], _n) for j, b in enumerate(RL.SOURCES)]
                                  for i, a in enumerate(RL.SOURCES)]))

    def _q(ids, sa=_sa, Ri=_Ri):
        # the GLS Q of the five ctrlonly estimates with the shared-case covariance: with S = D R D (D = diag(se)),
        # Q = u' R^-1 u - (w' R^-1 u)^2 / (w' R^-1 w), u = beta / se, w = 1 / se
        b = np.column_stack([sa[f'ctrlonly_{s}'].beta.reindex(ids).values for s in RL.SOURCES])
        se = np.column_stack([sa[f'ctrlonly_{s}'].se.reindex(ids).values for s in RL.SOURCES])
        ok = np.isfinite(b).all(1) & np.isfinite(se).all(1) & (se > 0).all(1)
        with np.errstate(all='ignore'):
            u, w = b / se, 1 / se
            qa = np.einsum('ij,jk,ik->i', u, Ri, u) - np.einsum('ij,jk,ik->i', w, Ri, u) ** 2 / \
                np.einsum('ij,jk,ik->i', w, Ri, w)
        return pd.Series(np.where(ok, qa, np.nan), index=ids)
    _h = _hq[(_hq.cohort == c) & (_hq.sources_used == len(RL.SOURCES))]
    check(len(_h) > 0 and np.allclose(_q(list(_h.variant_id)).values, _h.Q.values, rtol=1e-8, atol=1e-10),
          f'sources_null_q: {c}: the Q formula does not reproduce controls.heterogeneity.tsv')
    _p = _sa['pooled']
    _ids = [v for v in _nul.index if v in _p.index]
    _pp = _p.loc[_ids]
    _af = (_pp.af_case * _pp.n_case + _pp.af_ctrl * _pp.n_ctrl) / (_pp.n_case + _pp.n_ctrl)
    _d = pd.DataFrame({'Q': _q(_ids), 'maf_bin': np.minimum(_af, 1 - _af).map(PS.maf_bin).astype(str),
                       'is_indel': _nul.loc[_ids, 'is_indel'].values}, index=_ids).dropna()
    _d['p'] = SS.chi2.sf(_d.Q, len(RL.SOURCES) - 1)
    _cm = _d[_d.maf_bin.isin(_COMMON)]
    S['sources_null_q'][c] = {'n': int(len(_cm)), 'share_05': float((_cm.p < 0.05).mean()),
                              'share_thr': float((_cm.p < S['q_threshold']).mean())}
    S['sources_q_emp'][c] = {}
    for _g in S['q_markers']:
        _v = S['lead_marker'][_g]
        _r = _hq[(_hq.cohort == c) & (_hq.variant_id == _v)].iloc[0]
        _x = _p.loc[_v]
        _a = (_x.af_case * _x.n_case + _x.af_ctrl * _x.n_ctrl) / (_x.n_case + _x.n_ctrl)
        _ind = _v.startswith('chr') and (len(_v.split(':')[2]) != 1 or len(_v.split(':')[3]) != 1)
        _ref = _d[(_d.maf_bin == PS.maf_bin(min(_a, 1 - _a))) & (_d.is_indel == _ind)]
        check(abs(float(_r.p_Q) - S['sources_signal'][_g][c]['p_Q']) < 1e-12, f'sources_q_emp: {_g} {c}: Q P differs')
        S['sources_q_emp'][c][_g] = {'p': float(_r.p_Q), 'n_ref': int(len(_ref)),
                                     'p_emp': float((1 + (_ref.Q >= _r.Q).sum()) / (1 + len(_ref)))}
    for _sub in RL.PLATFORM_WITHOUT + RL.SOURCE_WITHOUT:
        _w = _sa[_sub].reindex(_ids)
        _sdc = np.sqrt(np.clip(1 - (_pp.se / _w.se) ** 2, 0, None))
        _u = (_w.z - _pp.beta / _w.se) / _sdc.where(_sdc > 0)
        _s = float(_u[np.isfinite(_u)].std())
        _sd['platforms' if _sub in RL.PLATFORM_WITHOUT else 'sources'].append(_s)
        if _sub in RL.SOURCE_WITHOUT:
            _und.append((float((~(_sdc > 0)).mean()), _sub.replace('ctrlwithout_', ''), c))
        if _sub == 'ctrlwithout_nagahama':
            _sd['nagahama'].append(_s)
S['loo_null_sd'] = {k: [min(v), max(v)] for k, v in _sd.items()}
_um = max(_und)
S['loo_null_undefined'] = {'share': _um[0], 'source': _um[1], 'cohort': _um[2]}
check(0.9 <= S['loo_null_sd']['platforms'][0] and S['loo_null_sd']['platforms'][1] <= 1.1,
      'loo_null_sd: the leave-one-out model no longer holds for the platforms; the notes say it does')
check(min(v['share_05'] for v in S['sources_null_q'].values()) > 0.05,
      'sources_null_q: the analytic source-Q P is no longer lenient; the slides say it is')

# ── reads at the two flagged leads, straight from the CRAMs (user 2026-11-06: "在27和28，我需要你查看cram原文件来
# 判明问题") ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# deck/cram_check/cram_scan.tsv (cram_check/README.md): every full-mainland genome, read only, joined to the pre-QC
# calls at the two leads. ADARB2 (chr10:1747167 CGTG>C, in a minisatellite of about 390 bp, 31-bp units; the carrier
# haplotype is one unit shorter, called as a 3-bp and a 28-bp deletion): reads of the minisatellite holding the carrier
# signature (16-mer GATCGCATGATGTGTA, absent from the reference window), whatever their alignment; a read-level carrier
# holds at least CRAM_MIN of them. LINC01146 (chr14:88019710 G>A): G and A reads with MAPQ and base quality >= 20; a
# read genotype from the A fraction (< 0.2, 0.2-0.8, > 0.8) where at least 8 such reads cover the site.
CR = pd.read_csv(os.path.join(HERE, 'cram_check', 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str})
check(len(CR) == S['n_case']['full_mainland'] + S['n_control']['full_mainland'], 'cram_scan: not every genome was read')
check(CR.error.isna().all(), 'cram_scan: a CRAM could not be read')
CRAM_MIN, LI_DEPTH = 4, 8
CR['ad_called'] = CR.ad_call.isin(['0/1', '1/1'])
CR['ad_read'] = CR.ad_alt >= CRAM_MIN
CR['li_depth'] = CR.li_ref_q + CR.li_alt_q
CR['li_frac'] = CR.li_alt_q / CR.li_depth.where(CR.li_depth > 0)
CR['li_read_gt'] = np.select([CR.li_frac < 0.2, CR.li_frac <= 0.8, CR.li_frac > 0.8], ['0/0', '0/1', '1/1'], 'NA')
_CG = {'controls': lambda d: d.status == 'control', 'T7': lambda d: d.platform == 'DNBSeq-T7 30x',
       'NovaSeq': lambda d: d.platform == 'NovaSeq 30x', 'G400RS 30x': lambda d: d.platform == 'DNBSeq-G400RS 30x',
       'G400RS 15x': lambda d: d.platform == 'DNBseq-G400RS 15x',
       'other cases': lambda d: (d.status == 'case') & (d.platform != 'DNBSeq-T7 30x'),
       'cases': lambda d: d.status == 'case'}
_DOS = {'0/0': 0, '0/1': 1, '1/1': 2}
S['cram'] = {'min_reads': CRAM_MIN, 'li_min_depth': LI_DEPTH, 'signature': 'GATCGCATGATGTGTA', 'genomes': int(len(CR)),
             'sig_bp': 16, 'unit_bp': 31, 'vntr': [1746881, 1747271], 'vntr_bp_about': 390, 'vntr_units_about': 13,
             'ref_units': 3, 'li_cuts': [0.2, 0.8], 'adarb2': {}, 'linc': {}}
check(S['cram']['unit_bp'] == 3 + S['adarb2']['full_mainland']['del28_bp'],
      'cram: the two called deletions no longer add up to one minisatellite unit')
for c in CK:
    x = CR[CR[f'in_{c}'] == 1]
    S['cram']['adarb2'][c], S['cram']['linc'][c] = {}, {}
    for g, f in _CG.items():
        y = x[f(x)]
        car = y[y.ad_called]
        S['cram']['adarb2'][c][g] = {
            'n': int(len(y)), 'called': int(len(car)), 'read': int(y.ad_read.sum()), 'called_read': int(car.ad_read.sum()),
            'missed': int((~y.ad_called & y.ad_read).sum()), 'called_pct': float(len(car) / len(y)),
            'read_pct': float(y.ad_read.mean()), 'zero': int((y.ad_alt == 0).sum()),
            'one_two': int(y.ad_alt.between(1, 2).sum()), 'three': int((y.ad_alt == 3).sum()),
            'sig_median_read': float(y[y.ad_read].ad_alt.median()) if y.ad_read.any() else None,
            'sig_range_read': [int(y[y.ad_read].ad_alt.min()), int(y[y.ad_read].ad_alt.max())] if y.ad_read.any() else None,
            'sig_q20_share': float(y[y.ad_read].ad_alt_q20.sum() / y[y.ad_read].ad_alt.sum()) if y.ad_read.any() else None,
            'rlen': float(y.ad_rlen_sum.sum() / y.ad_reads.sum()), 'reads_median': float(y.ad_reads.median()),
            'read_pct_k3': float((y.ad_alt >= 3).mean()),
            # review round 8: the called non-carriers alone (A-01: the controls' "3 signature reads" count included the
            # one called carrier below CRAM_MIN), and the soft-clipped share of the read carriers' signature reads (B-05)
            'nc_zero': int((~y.ad_called & (y.ad_alt == 0)).sum()),
            'nc_one_two': int((~y.ad_called & y.ad_alt.between(1, 2)).sum()),
            'nc_three': int((~y.ad_called & (y.ad_alt == 3)).sum()),
            'called_below_min': int((car.ad_alt < CRAM_MIN).sum()),
            'sig_clip_share': float(y[y.ad_read].ad_alt_clip.sum() / y[y.ad_read].ad_alt.sum()) if y.ad_read.any() else None}
        ok = y[y.li_depth >= LI_DEPTH]
        het, hr, ha = ok[ok.li_call == '0/1'], ok[ok.li_call == '0/0'], ok[ok.li_call == '1/1']
        S['cram']['linc'][c][g] = {
            'n': int(len(y)), 'depth_ok': int(len(ok)), 'concordant': int((ok.li_read_gt == ok.li_call).sum()),
            'concordance': float((ok.li_read_gt == ok.li_call).mean()),
            'het_ab_median': float(het.li_frac.median()) if len(het) else None,
            'het_ab_iqr': [float(het.li_frac.quantile(0.25)), float(het.li_frac.quantile(0.75))] if len(het) else None,
            'het_alt_fwd': float(het.li_alt_fwd_q.sum() / het.li_alt_q.sum()) if het.li_alt_q.sum() else None,
            'homref_alt': float(hr.li_alt_q.sum() / hr.li_depth.sum()) if len(hr) else None,
            'homalt_ref': float(ha.li_ref_q.sum() / ha.li_depth.sum()) if len(ha) else None,
            'depth_median': float(ok.li_depth.median()),
            'af_called': float(y.li_call.map(_DOS).sum() / (2 * len(y))),
            'af_reads': float(ok.li_read_gt.map(_DOS).sum() / (2 * len(ok)))}
_A = S['cram']['adarb2']
# review round 8 (A-05): the 150-bp platforms' calls among their read carriers
for c in CK:
    _p150 = ['controls', 'NovaSeq', 'G400RS 30x', 'G400RS 15x']
    _A[c]['p150'] = {'called_read': sum(_A[c][g]['called_read'] for g in _p150), 'read': sum(_A[c][g]['read'] for g in _p150)}
check(all(_A[c]['controls']['nc_zero'] + _A[c]['controls']['nc_one_two'] + _A[c]['controls']['nc_three'] +
          _A[c]['controls']['missed'] == _A[c]['controls']['n'] - _A[c]['controls']['called'] for c in CK),
      'cram: the controls called non-carriers do not add up')
check(all(_A[c]['T7']['called'] == 0 for c in CK), 'cram: a T7 case is now called an ADARB2 carrier')
check(all(_A[c]['controls']['called_read'] >= 0.99 * _A[c]['controls']['called'] for c in CK),
      'cram: the signature misses called ADARB2 carriers among the controls')
check(all(_A[c]['T7']['three'] == 0 and _A[c]['T7']['read_pct_k3'] == _A[c]['T7']['read_pct'] for c in CK),
      'cram: a DNBSeq-T7 count sits at the threshold; the slide says the threshold does not matter for T7')
check(all(S['cram']['adarb2'][c]['other cases']['called_read'] == S['cram']['adarb2'][c]['other cases']['called'] for c in CK),
      'cram: a called ADARB2 carrier among the other cases lacks the signature')
# where the ADARB2 deletion is lost (user 2026-11-06: "我们有没有办法确认是calling的问题？", then "26～29 重构"),
# cram_check/recall_check.py: every genome's alignment record (@PG), and the genotype at the lead after Q1 re-called the
# deep cases' reads at full depth (baseline arm, f = 1; full mainland, where every T7 case is)
_AL = pd.read_csv(os.path.join(HERE, 'cram_check', 'aligner.tsv'), sep='\t', dtype={'sample_id': str}, keep_default_na=False)
check(set(_AL.sample_id) == set(CR.sample_id) and (_AL.error == '').all() and (_AL.n_bwa.astype(int) == 1).all(),
      'aligner: not one bwa record for every scanned genome')
_cmd = _AL.groupby(['bwa_vn', 'bwa_cmd', 'samblaster_vn']).size()
S['cram']['aligner'] = {'genomes': int(len(_AL)), 'commands': int(len(_cmd)), 'platforms': int(_AL.platform.nunique()),
                        'bwa_vn': str(_cmd.index[0][0]), 'bwa_cmd': str(_cmd.index[0][1]),
                        'samblaster_vn': str(_cmd.index[0][2])}
check(S['cram']['aligner']['commands'] == 1, 'aligner: the genomes were aligned by more than one command')
_RC = pd.read_csv(os.path.join(HERE, 'cram_check', 'recall_adarb2.tsv'), sep='\t', dtype={'sample_id': str}, keep_default_na=False)
_RC['read'] = _RC.ad_alt >= CRAM_MIN
_RCC = ['0/1', '1/1']
_both = _RC[_RC.gt_recall != './.']
_t7r = _RC[(_RC.platform == 'DNBSeq-T7 30x') & _RC.read]
_o150 = _RC[(_RC.platform != 'DNBSeq-T7 30x') & _RC.read]
S['cram']['recall'] = {
    'genomes': int(len(_RC)), 't7': int((_RC.platform == 'DNBSeq-T7 30x').sum()),
    'other': int((_RC.platform != 'DNBSeq-T7 30x').sum()), 'both_called': int(len(_both)),
    'same': int((_both.gt_recall == _both.ad_call).sum()), 'missing': int((_RC.gt_recall == './.').sum()),
    't7_read': int(len(_t7r)), 't7_called': int(_t7r.gt_recall.isin(_RCC).sum()),
    't7_homref': int((_t7r.gt_recall == '0/0').sum()), 't7_missing': int((_t7r.gt_recall == './.').sum()),
    'other_read': int(len(_o150)), 'other_called_orig': int(_o150.ad_call.isin(_RCC).sum()),
    'other_called': int(_o150.gt_recall.isin(_RCC).sum())}
_RS = S['cram']['recall']
check(_RS['same'] == _RS['both_called'], 'recall: a re-called genotype at the ADARB2 lead differs from the original call')
check(_RS['t7_read'] == _A['full_mainland']['T7']['read'] and _RS['t7'] == _A['full_mainland']['T7']['n'],
      'recall: not every DNBSeq-T7 case was re-called')
check(_RS['t7_called'] == 0 and _RS['other_called'] > 0, 'recall: the re-call no longer misses the T7 read carriers only')
# HaplotypeCaller diagnostics at the lead (user 2026-11-06: "我们先看看HaplotypeCaller 的诊断文件"; then "26～29页重构
# … 说明其被flag的来源", 2026-11-07): cram_check/hc_debug (run_hc_debug.sh, parse_hc_log.py, hc_tables.py). 14 genomes
# re-called with the cohort's HaplotypeCaller command and output-only diagnostics: whether a candidate haplotype carries
# the one-unit deletion, the k-mer size of the assembly, the genotypes written; then the cohort's own records in the
# window before and after its site QC, and the G>A at 1747168 among the read carriers of every full-mainland genome
_HD = os.path.join(HERE, 'cram_check', 'hc_debug')
_HG = pd.read_csv(os.path.join(_HD, 'hc_genotypes.tsv'), sep='\t', dtype={'sample_id': str, 'kmer_used': str}, keep_default_na=False)
_HR = pd.read_csv(os.path.join(_HD, 'cohort_records.tsv'), sep='\t').set_index('pos')
_HC = pd.read_csv(os.path.join(_HD, 'cohort_ga.tsv'), sep='\t').set_index('group')
_HS = pd.read_csv(os.path.join(_HD, 'hc_summary.tsv'), sep='\t', dtype={'sample_id': str})
_role = {'t7_read': _HG.role == 'T7 read carrier', 'p150_called': _HG.role.str.endswith('called carrier'),
         'p150_missed': _HG.role == 'HiSeq X missed read carrier', 't7_non': _HG.role == 'T7 non-carrier'}
check(sum(int(m.sum()) for m in _role.values()) == len(_HG), 'hc_debug: a diagnosed genome has no role')
_ga_het = lambda r: [int(r[f'{p}_gq']) for p in ('ga_1747168', 'ga_1747199') if r[f'{p}_gt'] == '0/1']  # noqa: E731
_hcg = {}
for k, m in _role.items():
    g = _HG[m]
    gq = [max(_ga_het(r)) for _, r in g.iterrows() if _ga_het(r)]
    _hcg[k] = {'n': int(len(g)), 'del_hap': int((g.deletion_haplotypes > 0).sum()),
               'lead_called': int(g.lead_gt.isin(['0/1', '1/1']).sum()), 'ga_het': len(gq),
               'ga_gq_min': min(gq) if gq else None,
               'reads_median': float(_HS[_HS.sample_id.isin(g.sample_id)].reads.median())}
_k = sorted({int(x) for v in _HS.kmers_used.astype(str) for x in v.split(',')})
_kd = sorted({int(x.split(':')[0]) for v in _HS.kmers_dropped for x in v.split()})
check(_k == [85] and all('non-unique' in x for v in _HS.kmers_dropped for x in v.split()),
      'hc_debug: not every assembly over the lead used k = 85 after dropping smaller k for non-unique k-mers')
_cr = lambda p: _HR.loc[p]  # noqa: E731
S['cram']['hc'] = {
    'genomes': int(len(_HG)), 'carriers': int((~_role['t7_non']).sum()), 'ga_pos': [1747168, 1747199], 'k_used': _k[0], 'k_dropped': _kd, 'kmers_per_read': {'100': 100 - _k[0] + 1, '150': 150 - _k[0] + 1},
    'groups': _hcg,
    'cohort': {'vqslod_min': 10, 'mq_min': 58.75,
               'del_vqslod': float(_cr(1747167).vqslod), 'del_kept': int(_cr(1747167).kept_after_site_qc),
               'del28_vqslod': float(_cr(1747171).vqslod),
               'ga_vqslod': float(_cr(1747168).vqslod), 'ga_mq': float(_cr(1747168).mq), 'ga_ac': int(_cr(1747168).ac),
               'ga_kept': int(_cr(1747168).kept_after_site_qc), 'ga1747199_in_calls': int(1747199 in _HR.index),
               't7_missed': int(_HC.loc['DNBSeq-T7 cases', 'read_carriers_missed']),
               't7_missed_ga': int(_HC.loc['DNBSeq-T7 cases', 'missed_with_ga']),
               'ctrl_missed': int(_HC.loc['HiSeq X controls', 'read_carriers_missed']),
               'ctrl_missed_ga': int(_HC.loc['HiSeq X controls', 'missed_with_ga']),
               'g15_missed': int(_HC.loc['G400RS 15x cases', 'read_carriers_missed']),
               'g15_missed_ga': int(_HC.loc['G400RS 15x cases', 'missed_with_ga']),
               'missed': int(_HC.loc['all', 'read_carriers_missed']), 'missed_ga': int(_HC.loc['all', 'missed_with_ga']),
               'called': int(_HC.loc['all', 'called_carriers']), 'called_ga': int(_HC.loc['all', 'called_with_ga']),
               'noncarrier_ga': int(_HC.loc['all', 'noncarriers_with_ga']),
               # review round 9 (B-06): AC counts every genome of the cohort's VCF; the carriers above are full mainland
               'cohort_genomes': int(_cr(1747168).an) // 2, 'ga_carriers_fm': int(_HC.loc['all', 'ga_carriers'])}}
# review round 9 (A-20): the signature reads of the diagnosed genomes, by group
for _k, _m in _role.items():
    S['cram']['hc']['groups'][_k]['sig_range'] = [int(_HG[_m].signature_reads.min()), int(_HG[_m].signature_reads.max())]
_H, _HCo = S['cram']['hc']['groups'], S['cram']['hc']['cohort']
check(_H['t7_read']['del_hap'] == 0 and _H['p150_called']['del_hap'] == _H['p150_called']['n'] == _H['p150_called']['lead_called']
      and _H['p150_missed']['del_hap'] == 0 and _H['t7_read']['ga_het'] == _H['t7_read']['n'],
      'hc_debug: the deletion haplotype is no longer assembled in the called 150-bp carriers only')
check(_HCo['ga_kept'] == 0 and _HCo['del_kept'] == 1 and _HCo['ga_vqslod'] < _HCo['vqslod_min'] < _HCo['del_vqslod']
      and _HCo['called_ga'] == 0, 'hc_debug: the G>A record no longer fails the site QC that the deletion passes')
_AF = S['cram']['adarb2']['full_mainland']
check(_HCo['t7_missed'] == _AF['T7']['read'] and _HCo['ctrl_missed'] == _AF['controls']['missed'] and
      _HCo['called'] == _AF['cases']['called'] + _AF['controls']['called'],
      'hc_debug: the cohort counts disagree with cram_scan')
# HaplotypeCaller census at the lead (user 2026-11-07: "我们有没有什么办法取样取的更加系统有说服力"): cram_check/hc_census
# (make_samples.py, run_hc_census.sh, census_tables.py): the pilot's command on every full-mainland genome that carries
# the deletion in its reads or in the calls, and on 20 random non-carriers per platform (fixed seed). Counts only here:
# the per-sample table (hc_census.tsv, sample IDs) is internal and git-ignored; stats.json may be published.
_CD = os.path.join(HERE, 'cram_check', 'hc_census')
_CT = pd.read_csv(os.path.join(_CD, 'hc_census.tsv'), sep='\t', dtype={'sample_id': str, 'kmers_used': str},
                  keep_default_na=False)
_CGRP = {'called': ['HiSeq X controls, called', '150-bp cases, called'],
         'called_ctrl': ['HiSeq X controls, called'], 'called_case': ['150-bp cases, called'],
         'missed': ['HiSeq X controls, missed', '150-bp cases, missed'],
         'missed_ctrl': ['HiSeq X controls, missed'], 'missed_case': ['150-bp cases, missed'],
         't7': ['DNBSeq-T7 read carriers'], 'nosig': ['Called, no signature read'],
         'non': ['Non-carriers, 150-bp platforms', 'Non-carriers, DNBSeq-T7'],
         'non_150': ['Non-carriers, 150-bp platforms'], 'non_t7': ['Non-carriers, DNBSeq-T7']}


def _cgroup(gs):
    g = _CT[_CT.group.isin(gs)]
    return {'n': int(len(g)), 'del': int(g.deletion_assembled.sum()), 'lead': int(g.rerun_lead_carrier.sum()),
            'ga': int(g.rerun_ga_carrier.sum()), 'repro': int(g.reproduces_cohort.sum()),
            'active': int(g.lead_region_active.sum()),
            'sig_median': float(g.signature_reads.median()) if len(g) else None,
            'sig_range': [int(g.signature_reads.min()), int(g.signature_reads.max())] if len(g) else None,
            'reads_median': float(g.region_reads.median()) if len(g) else None}


_ck = sorted({int(x) for v in _CT.kmers_used for x in str(v).split(';') if x})
S['cram']['hc_census'] = {
    'genomes': int(len(_CT)), 'carriers': int((_CT.role != 'non-carrier, random').sum()),
    'noncarriers_per_platform': 20, 'seed': 20261107, 'k_used': _ck,
    'groups': {k: _cgroup(v) for k, v in _CGRP.items()},
    'by_platform': {p: _cgroup_p for p, _cgroup_p in
                    ((p, {'n': int(len(g)), 'del': int(g.deletion_assembled.sum()),
                          'called': int(g.role.eq('read carrier, called').sum() + g.role.eq('called, no signature').sum()),
                          'missed': int(g.role.eq('read carrier, missed').sum())})
                     for p, g in _CT[_CT.role != 'non-carrier, random'].groupby('platform'))}}
_CS = S['cram']['hc_census']
# the census is every carrier of cram_scan and 20 non-carriers per platform, and nothing else
check(_CS['carriers'] == _AF['cases']['read'] + _AF['controls']['read'] + int((_CT.role == 'called, no signature').sum())
      and _CS['groups']['t7']['n'] == _AF['T7']['read'] and _CS['groups']['missed_ctrl']['n'] == _AF['controls']['missed']
      and _CS['groups']['called']['n'] + _CS['groups']['nosig']['n'] == _AF['cases']['called'] + _AF['controls']['called']
      and _CS['groups']['non']['n'] == 5 * _CS['noncarriers_per_platform'],
      'hc_census: the census genomes disagree with cram_scan')
check(sum(_CS['groups'][k]['n'] for k in ('called', 'missed', 't7', 'nosig', 'non')) == _CS['genomes'],
      'hc_census: a genome has no group')
# ToMMo 61KJPN at the signature (user 2026-11-07, round 12 item 10; B-12): the 16-bp signature arises from the one-unit
# deletion or from a G>A in a full-length unit; cram_check/tommo_signature.py (allele counts only)
_TS = pd.read_csv(os.path.join(HERE, 'cram_check', 'tommo_signature.tsv'), sep='\t').set_index('allele')
S['cram']['tommo_sig'] = {a: {'ac': int(r.ac), 'an': int(r.an), 'af': float(r.af), 'carriers_expected': float(r.carriers_expected),
                              'filter': str(r['filter'])} for a, r in _TS.iterrows()}
_tg = S['cram']['tommo_sig']
check(abs(_tg['deletion_3bp']['af'] - S['lead_af']['full_mainland']['ADARB2']['reference | ToMMo 61KJPN']['af']) < 0.002
      and _tg['ga_1747168']['af'] < 0.01 * _tg['deletion_3bp']['af'] * 10,
      "tommo_sig: the deletion no longer matches the lead's ToMMo frequency, or the G>A is no longer rare beside it")
# what the slide says: the re-run reproduces every call; the deletion is reconstructed in every called carrier and in no
# other genome (missed 150-bp carriers, T7 carriers, random non-carriers); k = 85 in every assembly over the lead
_CG = _CS['groups']
_CS['reproduced'] = int(_CT.reproduces_cohort.sum())
_CS['assembled_is_called'] = int(((_CT.deletion_assembled == 1) == _CT.cohort_call.astype(str).str.contains('[1-9]')).sum())
check(_CS['reproduced'] == _CS['genomes'] == _CS['assembled_is_called'] and _CS['k_used'] == [85]
      and _CG['called']['del'] == _CG['called']['n'] and _CG['nosig']['del'] == _CG['nosig']['n']
      and _CG['missed']['del'] == 0 and _CG['t7']['del'] == 0 and _CG['non']['del'] == 0,
      'hc_census: the reconstruction no longer matches the calls')

# ── LINC01146 and the G400RS 15x cases (user 2026-11-06: "LINC01146 应该是G400RS 15X的问题，你看看") ───────────────
# The flag (leaving out the DNBSeq-T7 cases) read against each case platform: the platform strata of the pipeline
# (09.platform_assoc/<cohort>/only_<platform>.saige.tsv; shared-control GLS as platform_stats.py), each stratum's
# departure from the pooled estimate, the strata combined without T7 and without both T7 and G400RS 15x, and Q with
# each platform left out; then the G400RS 15x genomes: independent fragments at the lead (cram_check/linc_fragments.tsv)
# and genome-wide heterozygosity and missingness on the common null SNVs (cram_check/het_by_sample.tsv).
from scipy.stats import mannwhitneyu as _MW  # noqa: E402
_LV = S['lead_marker']['LINC01146']
_KS = ['t7', 'novaseq', 'g400rs30', 'g400rs15']
S['linc2'] = {}
for c in CK:
    _cnt = pd.read_csv(RL.subset_counts_file(R, c, 'pooled'), sep='\t').set_index('subset')
    _n0 = int(_cnt.loc['pooled', 'controls_random'])
    _sa = RL.read_saige_subsets(R, c, ['pooled'] + [f'only_{k}' for k in _KS])
    _nk = {k: int(_cnt.loc[f'only_{k}', 'cases_random']) for k in _KS}
    _bp, _sp = float(_sa['pooled'].loc[_LV, 'beta']), float(_sa['pooled'].loc[_LV, 'se'])

    def _gls(ks, sa=_sa, nk=_nk, n0=_n0):
        b = np.array([sa[f'only_{k}'].loc[_LV, 'beta'] for k in ks])
        se = np.array([sa[f'only_{k}'].loc[_LV, 'se'] for k in ks])
        Rm = np.array([[1.0 if i == j else PS.rho_shared(nk[a], nk[bb], n0) for j, bb in enumerate(ks)]
                       for i, a in enumerate(ks)])
        mu, sem, q = PS.gls(b, Rm * np.outer(se, se))
        return float(mu), float(sem), float(q)
    L2 = {'pooled_or': math.exp(_bp), 'pooled_z': _bp / _sp, 'pooled_se': _sp, 'strata': {}, 'combos': {}, 'q': {}}
    for k in _KS:
        b, se = float(_sa[f'only_{k}'].loc[_LV, 'beta']), float(_sa[f'only_{k}'].loc[_LV, 'se'])
        L2['strata'][k] = {'n': _nk[k], 'or': math.exp(b), 'lo': math.exp(b - 1.959964 * se), 'hi': math.exp(b + 1.959964 * se),
                           'departure': (b - _bp) / se}
    # review round 8 (B-01): each stratum against the other cases, the pipeline's difference test (platform_stats: the
    # stratum against the fit without it, with the shared-control covariance), signed; all cases are 73% T7
    _lp = tsv('10.platform', c, 'platform.lopo.tsv', dtype={'variant_id': str})
    _lp = _lp[_lp.variant_id == _LV].set_index('left_out')
    for k in _KS:
        L2['strata'][k]['z_other'] = float(_lp.loc[k, 'diff'] / _lp.loc[k, 'se_diff'])
        L2['strata'][k]['p_other'] = float(_lp.loc[k, 'p_diff'])
    check(all(abs(L2['strata'][k]['p_other'] - S['platform']['LINC01146'][c][f'pdiff_{k}']) < 1e-9 for k in _KS),
          f'linc2: {c}: the strata difference tests differ from platform.signals')
    # review round 8 (A-03): what the association reads at the lead. The subsets' SAIGE allele frequencies count each
    # genotype set missing at QC as its best guess (one A copy here); the genotypes set missing per group
    _q = lpq[(lpq.cohort == c) & (lpq.signal == 'LINC01146')].set_index('group')
    L2['qc_missing'] = {g: {'n': int(_q.loc[grp, 'n']), 'missing': int(_q.loc[grp, 'set_missing_by_filters'])}
                        for g, grp in (('cases', 'case | all'), ('controls', 'control | all'),
                                       ('t7', 'case | DNBSeq-T7 30x'), ('g400rs15', 'case | DNBseq-G400RS 15x'))}
    L2['input_af'] = {'cases': float(_sa['pooled'].loc[_LV, 'af_case']), 'controls': float(_sa['pooled'].loc[_LV, 'af_ctrl']),
                      't7': float(_sa['only_t7'].loc[_LV, 'af_case']), 'g400rs15': float(_sa['only_g400rs15'].loc[_LV, 'af_case'])}
    for g, grp in (('controls', 'control | all'), ('t7', 'case | DNBSeq-T7 30x'), ('g400rs15', 'case | DNBseq-G400RS 15x')):
        _a, _m = S['lead_af'][c]['LINC01146'][grp], L2['qc_missing'][g]
        check(abs(L2['input_af'][g] - (round(2 * _a['n'] * _a['af']) + _m['missing']) / (2 * _m['n'])) < 1e-3,
              f'linc2: {c}: {g}: the input frequency is not the QC\'d genotypes plus one A copy per missing genotype')
    _cl = S['cram']['linc'][c]
    L2['af_diff'] = {'input': L2['input_af']['cases'] - L2['input_af']['controls'],
                     'reads': _cl['cases']['af_reads'] - _cl['controls']['af_reads'],
                     'pre_qc': _cl['cases']['af_called'] - _cl['controls']['af_called']}
    # review round 8, U5 (A-03; user 2026-11-06 "U1-U6 按照最合理的建议修改"): a descriptive reading with the genotypes
    # set missing at QC left out (followup.gwas_af.tsv counts the QC'd genotypes only), and what those genotypes carry
    # in the pre-QC calls: the CRAM table's pre-QC A count over every genome less the QC'd genotypes' count
    _ca = S['lead_af'][c]['LINC01146']
    L2['called_af'] = {'cases': float(_ca['case | all']['af']), 'controls': float(_ca['control | all']['af'])}
    L2['af_diff']['called'] = L2['called_af']['cases'] - L2['called_af']['controls']
    L2['missing_pre'] = {}
    for g, grp in (('cases', 'case | all'), ('controls', 'control | all')):
        _m = L2['qc_missing'][g]
        check(int(_ca[grp]['n']) + _m['missing'] == _m['n'] == _cl[g]['n'],
              f'linc2: {c}: {g}: the QC\'d and the missing genotypes do not add up to the group')
        _apre, _aqc = 2 * _m['n'] * _cl[g]['af_called'], 2 * _ca[grp]['n'] * _ca[grp]['af']
        check(abs(_apre - round(_apre)) < 1e-6 and abs(_aqc - round(_aqc)) < 1e-6,
              f'linc2: {c}: {g}: an A count is not whole (pre-QC calls over every genome, QC\'d genotypes)')
        L2['missing_pre'][g] = {'a': int(round(_apre) - round(_aqc)), 'alleles': 2 * _m['missing'],
                                'with_a': int(_q.loc[grp, 'set_missing_with_alt'])}
        check(0 <= L2['missing_pre'][g]['a'] <= L2['missing_pre'][g]['alleles'],
              f'linc2: {c}: {g}: the genotypes set missing carry an impossible A count')
    for key, ks in (('without_t7', ['novaseq', 'g400rs30', 'g400rs15']), ('without_t7_g400rs15', ['novaseq', 'g400rs30']),
                    ('without_g400rs15', ['t7', 'novaseq', 'g400rs30'])):
        mu, sem, _ = _gls(ks)
        sdc = math.sqrt(max(1e-12, 1 - (_sp / sem) ** 2))
        L2['combos'][key] = {'or': math.exp(mu), 'se': sem, 'z': mu / sem, 'z_exp': _bp / sem,
                             'chance': float(SS.norm.cdf((mu / sem - _bp / sem) / sdc))}
    for key, ks in (('all', _KS), *[(f'without_{k}', [x for x in _KS if x != k]) for k in _KS]):
        q = _gls(ks)[2]
        L2['q'][key] = {'q': q, 'p': float(SS.chi2.sf(q, len(ks) - 1))}
    # the GLS of the strata stands in for the pipeline's own fit without T7 (checked)
    _pl = loo[(loo.signal == 'LINC01146') & (loo.check == 'lopo') & (loo.cohort == c) & (loo.left_out == 't7')].iloc[0]
    check(abs(L2['combos']['without_t7']['z'] - abs(_pl.z_observed)) < 0.15,
          f'linc2: {c}: the strata combined without T7 drift from the pipeline\'s fit')
    S['linc2'][c] = L2
check(all(min(S['linc2'][c]['strata'].values(), key=lambda v: v['departure']) is S['linc2'][c]['strata']['g400rs15']
          and S['linc2'][c]['strata']['g400rs15']['or'] < 1 for c in CK),
      'linc2: G400RS 15x is no longer the stratum furthest below the pooled estimate')
# review round 8, U6 (B-08; user 2026-11-06 "U1-U6 按照最合理的建议修改"): the cases' recruiting institution, from the
# sample sheet (Providing_Institution), per cohort; counts only. The registration study is already split in
# followup.gwas_af.tsv (lead_af: '<platform> | old study' / '| new study')
_inf = pd.read_excel(RL.K.INFO_XLSX, usecols=['ID_JHRPv6', 'Providing_Institution'], dtype=str).dropna(subset=['ID_JHRPv6'])
check(not _inf.ID_JHRPv6.duplicated().any(), 'institution: the sample sheet repeats an ID')
_inf = _inf.set_index('ID_JHRPv6')
_cs = stt[stt.status == 'case'].copy()
_cs['inst'] = _cs.sample_id.map(_inf.Providing_Institution)
check(_cs.inst.notna().all(), 'institution: a case lacks Providing_Institution in the sample sheet')
_top = _cs[_cs.in_full_mainland == '1'].inst.value_counts().index[0]
S['institution'] = {}
for c in CK:
    _x = _cs[_cs[f'in_{c}'] == '1']
    _o = _x[_x.inst != _top]
    S['institution'][c] = {'cases': int(len(_x)), 'main': int(len(_x) - len(_o)), 'others': int(len(_o)),
                           'other_institutions': int(_o.inst.nunique()), 'other_platforms': sorted(set(_o.platform))}
    check(S['institution'][c]['cases'] == S['linc2'][c]['qc_missing']['cases']['n'],
          f'institution: {c}: the sample table\'s cases differ from the association\'s')
check(all(set(S['institution'][c]['other_platforms']) <= {'NovaSeq 30x'} for c in CK),
      'institution: a case from another institution was sequenced on a DNBSeq platform; the notes say every DNBSeq '
      'case comes from one institution')
_FR = pd.read_csv(os.path.join(HERE, 'cram_check', 'linc_fragments.tsv'), sep='\t', dtype={'sample_id': str})
_FR = _FR.merge(CR[['sample_id', 'li_call', 'li_dp']], on='sample_id')
_g15 = _FR[_FR.platform == 'DNBseq-G400RS 15x']
_hom = _g15[_g15.li_call.isin(['0/0', '1/1'])]
S['g15'] = {'n': int(len(_g15)), 'reads_median': float(_g15.reads.median()), 'dup': int(_g15.dup_flagged.sum()),
            'reads': int(_g15.reads.sum()), 'hom': int(len(_hom)), 'hom_frag_min': int(_hom.fragments.min()),
            'hom_other_allele': int(((_hom.li_call == '0/0') & (_hom.A > 0)).sum() + ((_hom.li_call == '1/1') & (_hom.G > 0)).sum()),
            'calls': {k: int((_g15.li_call == k).sum()) for k in ('0/0', '0/1', '1/1')}}
_n = S['g15']['calls']
_p = (_n['0/1'] + 2 * _n['1/1']) / (2 * S['g15']['n'])
S['g15']['het_expected'] = float(2 * _p * (1 - _p) * S['g15']['n'])
S['g15']['het_dropout_p_min'] = float(2 * 0.5 ** S['g15']['hom_frag_min'])
# review round 8 (A-03, B-07): the genotypes the association reads. The cohort's QC sets a genotype missing below
# DP 8 (const dp); the two thinnest G400RS 15x genomes fall there and enter the test as their best guess
_QC_DP = 8
_g15q = _g15[_g15.li_dp >= _QC_DP]
_homq = _g15q[_g15q.li_call.isin(['0/0', '1/1'])]
S['g15']['qc'] = {'n': int(len(_g15q)), 'missing': int(len(_g15) - len(_g15q)),
                  'missing_dp': sorted(int(x) for x in _g15[_g15.li_dp < _QC_DP].li_dp),
                  'calls': {k: int((_g15q.li_call == k).sum()) for k in ('0/0', '0/1', '1/1')},
                  'frag_min': int(_g15q.fragments.min()), 'hom_frag_min': int(_homq.fragments.min())}
_nq = S['g15']['qc']['calls']
_pq = (_nq['0/1'] + 2 * _nq['1/1']) / (2 * S['g15']['qc']['n'])
S['g15']['qc']['het_expected'] = float(2 * _pq * (1 - _pq) * S['g15']['qc']['n'])
S['g15']['qc']['het_dropout_p_max'] = float(2 * 0.5 ** S['g15']['qc']['hom_frag_min'])
check(S['g15']['qc']['missing'] == int(lpq[(lpq.cohort == 'full_mainland') & (lpq.signal == 'LINC01146') &
                                           (lpq.group == 'case | DNBseq-G400RS 15x')].set_missing_by_filters.iloc[0]),
      'g15: the genomes below DP 8 are not the genotypes QC set missing')
_HB = pd.read_csv(os.path.join(HERE, 'cram_check', 'het_by_sample.tsv'), sep='\t', dtype={'sample_id': str})
_HB['g'] = np.where(_HB.status == 'control', 'controls', _HB.platform)
_ref = _HB[_HB.g == 'controls']
S['het_genome'] = {g: {'n': int(len(x)), 'F': float(x.F.median()), 'missing': float(x.F_MISS.median()),
                       'depth': float(x.observed_depth.median()),
                       'p_F': None if g == 'controls' else float(_MW(x.F, _ref.F).pvalue),
                       'p_missing': None if g == 'controls' else float(_MW(x.F_MISS, _ref.F_MISS).pvalue)}
                   for g, x in _HB.groupby('g')}
_lo = _ref[_ref.observed_depth.between(14, 16)]
S['het_genome']['controls_14_16x'] = {'n': int(len(_lo)), 'F': float(_lo.F.median()), 'missing': float(_lo.F_MISS.median())}
S['het_genome_snvs'] = 6000

# ── constants of the design (docs/METHODS.md, evidence.rules.tsv) ──────────
# the evidence-rule thresholds come from the registry (scripts/checks.py); the rest are the pipeline's settings
RK = DF.CK
S['const'] = {'z_exp_min': RK.LOO_Z_EXPECTED, 'z_obs_max': RK.LOO_Z_OBSERVED, 'depth_drop': RK.DEPTH_DZ,
              'cohorts_to_flag': RK.MIN_COHORTS, 'depth_cohorts_to_flag': RK.min_cohorts('depth'), 'cohorts': len(DF.COHORTS),
              'maf_q': RK.Q_MIN_MAF, 'gw': 5e-8, 'sg': 1e-5, 'pcs': 10,
              'gwas_variants': 168, 'loci': 16, 'alt_leads': 3, 'cs_members': 149, 'null_snps': 20000,
              'proxy_r2': 0.8, 'array_dz': 3, 'array_p': 0.05,
              'imputed_r2': 0.8, 'g6_conc': 0.995, 'g6_callrate_pp': 0.5, 'g8_lo': 0.95, 'g8_hi': 1.05,
              'g8_any_lo': 0.85, 'g8_any_hi': 1.15, 'g5_tol': 0.01, 'g10_tol': 1e-6, 'interval_pad_kb': 10,
              'shard_bp': 5000, 'gq': 20, 'dp': 8, 'ab_lo': 0.2, 'ab_hi': 0.8, 'hlahd': '1.7.1', 'gatk_call': '4.1.4.0',
              'gatk_qc': '4.6.2.0', 'power_80': 80, 'one_in_five': 5, 'mq': 20, 'bq': 10, 'readlen_reads': 1000,
              'probe_windows': 10, 'probe_mb': 1, 'deck_date': '2026-11-07', 'alt_reads_min': 2,
              # a leave-one-out reading with |z_exp| below this has almost no information left (named in the figure)
              'void_z': 0.5,
              # the testability bar is 1.96 plus this margin: an unchanged effect at the bar reaches 1.96 with 80% power
              'z_margin': round(RK.LOO_Z_EXPECTED - RK.LOO_Z_OBSERVED, 2)}

# ── references (verified: Crossref and PubMed records, 2026-10-31) ─────────
S['refs'] = json.load(open(os.path.join(HERE, 'references.json')))

# ── the design section: the threat model, coverage and registry of scripts/design_figures.py ────────────────────
S['design_figs'] = json.load(open(os.path.join(HERE, '..', 'docs', 'design', 'figures.json')))
S['design_threats'] = DF.threats(N)
S['design_coverage'] = DF.COVERAGE
S['design_experiments'] = [dict(zip(('key', 'name', 'checks', 'verb'), e)) for e in DF.EXPERIMENTS]
S['design_rows'] = [dict(zip(('label', 'what', 'roles', 'note_cases', 'note_controls'), r)) for r in DF.coverage_rows(N)]
S['design_roles'] = {k: {'fill': v[0], 'text': v[1], 'label': v[2]} for k, v in DF.ROLE_STYLE.items()}
S['registry'] = {'checks': [{'key': c.key, 'title': c.title, 'artifact': c.artifact, 'samples': c.samples,
                             'readout': c.readout, 'rule': c.rule, 'tables': c.tables, 'methods': c.methods}
                            for c in DF.CK.CHECKS],
                 'gates': [{'key': g.key, 'what': g.what, 'rule': g.rule, 'process': g.process} for g in DF.CK.GATES],
                 'processes': [{'name': p.name, 'stage': p.stage, 'feeds': p.feeds} for p in DF.PROCESSES],
                 'readings': DF.CK.READINGS}
S['warnings'] = WARN

# sample IDs never enter stats.json (user 2026-11-07, round 12 item 7: the repository is public, and git would commit a
# .json file). The DQB1:57:V trace keeps its cases under pseudonyms ("traced case 1", … in the order of the sorted IDs):
# the deck joins on them and names a traced case by its platform only. Then every sample ID of the sample table is
# searched for in the serialized file.
_vid = sorted(set(S['v57_traced']) | set(S['v57_changed']) | set(S['v57_ambiguous']))
_ps = {s: f'traced case {i + 1}' for i, s in enumerate(_vid)}
S['v57_traced'] = [_ps[s] for s in S['v57_traced']]
S['v57_changed'] = [_ps[s] for s in S['v57_changed']]
S['v57_ambiguous'] = [_ps[s] for s in S['v57_ambiguous']]
S['v57_traced_platform'] = {_ps[s]: p for s, p in S['v57_traced_platform'].items()}
for _r in S['v57_trace']:
    _r['sample_id'] = _ps[_r['sample_id']]
_ids = set(tsv('00.samples', 'sample_table.tsv', dtype=str).sample_id) | set(
    pd.read_csv(os.path.join(HERE, 'cram_check', 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str}).sample_id)
_js = json.dumps(clean(S), ensure_ascii=False)
_hit = {m.group(0) for m in re.finditer(r'(?<![\w])(?:' + '|'.join(re.escape(x) for x in sorted(_ids, key=len, reverse=True))
                                        + r')(?![\w])', _js)}
check(not _hit, f'stats.json would hold {len(_hit)} sample ID(s)')
json.dump(clean(S), open(os.path.join(HERE, 'stats.json'), 'w'), indent=1, ensure_ascii=False)
print(f'wrote stats.json ({S["roots"]["mode"]}):', len(json.dumps(clean(S))), 'bytes;', S['n_signals'], 'signals,',
      S['flag_cells'], 'flag cells; flagged:', S['flag_signals'])
for w in WARN:
    print('  warning:', w)
