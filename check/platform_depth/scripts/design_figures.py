#!/usr/bin/env python3
"""The design figures of check/platform_depth (documentation only: nothing here is staged by the workflow, and
nothing under results/ is written). Every number is read from the registry (scripts/checks.py), the published
tables (results/), the run traces (results/_run_info/*/trace.txt) and nextflow.config; none is typed here.

  numbers()        the numbers the figures print; deck/collect_stats.py imports it, so the deck's design section
                   and these figures print the same values
  workflow_edges() the process graph: the channels of platform_depth.nf (cached stages and the report stage's own
                   channels) and the published folders the report stage reads (the registry's `reads`)
  figures          docs/design/<name>.png (300 dpi) with their sources (<name>.dot for graphviz):
                     threat_model     what each artifact would fake, the prediction it makes, the experiment that
                                      tests it (matplotlib)
                     master_workflow  data -> fidelity gates -> four experiments -> readouts -> evidence matrix ->
                                      decision, with the key numbers (graphviz)
                     numbers_flow     the units at every step, CONSORT style, and the compute used (graphviz)
                     coverage         which genomes each experiment changes, to scale (matplotlib)
                     process_dag      the documented processes by stage, gates marked, published folders (graphviz)
                   docs/design/figures.json: natural size and smallest font of each figure (the deck places them
                   at a size that keeps every text >= 9 pt)
  --check-doc      docs/DESIGN.md must be what this script writes from the current numbers (exit 1 otherwise)

Where the numbers come from (environment, read once at import):
  PD_RESULTS   the results root (default: <component>/results)
  PD_TABLES    optional: one flat folder holding the report stage's tables in their final format (ds.*, array.*,
               controls.*, followup.*, evidence.*), read instead of their published folders under PD_RESULTS (a
               development copy of the report stage, e.g. .dev_report/tabs); every other table comes from PD_RESULTS

Usage (from the component folder):
  PYTHONDONTWRITEBYTECODE=1 /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -B scripts/design_figures.py
  PYTHONDONTWRITEBYTECODE=1 /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -B scripts/design_figures.py --check-doc
graphviz: /usr/bin/dot; fonts: Liberation Sans through deck/fonts.conf (the deck's PDF font for DengXian).
"""
import argparse
import glob
import json
import math
import os
import re
import subprocess
import sys
import textwrap
from collections import OrderedDict, defaultdict
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = Path(os.environ.get('PD_RESULTS') or ROOT / 'results').resolve()
TABLES = Path(os.environ['PD_TABLES']).resolve() if os.environ.get('PD_TABLES') else None
OUT = ROOT / 'docs' / 'design'
sys.path.insert(0, str(HERE))
import checks as CK  # noqa: E402  the registry (data only)

# the published folder of each report table, by its file prefix (docs/OUTPUTS.md)
REPORT_FOLDER = {'ds': ('07.downsampling', 'downsampling', 'summary'), 'array': ('11.array',),
                 'controls': ('12.controls',), 'followup': ('13.followup',), 'evidence': ('14.evidence',)}
# processes that run in the workflow but are not part of the documented design: left out of the process graph, the
# process counts and the compute (their tables are never read here)
UNDOCUMENTED = {'HLA_REF_SHIFT'}
PROCESSES = [p for p in CK.PROCESSES if p.name not in UNDOCUMENTED]


def report_file(name):
    """a report table: <PD_TABLES>/<name> when PD_TABLES is set, else its published folder under the results root"""
    if TABLES is not None:
        return TABLES / name
    return RES.joinpath(*REPORT_FOLDER[name.split('.')[0]], name)

COHORTS = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
CSHORT = {'narrow_mainland': 'N', 'intermediate_mainland': 'I', 'full_mainland': 'F'}
CASE_PLATFORMS = ['DNBSeq-T7 30x', 'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x']
PLAT_SHORT = {'DNBSeq-T7 30x': 'DNBSeq-T7', 'NovaSeq 30x': 'NovaSeq', 'DNBSeq-G400RS 30x': 'G400RS 30x',
              'DNBseq-G400RS 15x': 'G400RS 15x', 'HiSeqX 15x': 'HiSeq X'}
CONTROL_PLATFORM = 'HiSeqX 15x'
SOURCES = ['nagahama', 'bbj', 'htlv1', 'acc', 'lungca']
SRC_LABEL = {'nagahama': 'Nagahama', 'bbj': 'BBJ', 'htlv1': 'HTLV-1', 'acc': 'ACC', 'lungca': 'Lung cancer'}
FONTS_CONF = ROOT / 'deck' / 'fonts.conf'
DOT = '/usr/bin/dot'
FONT = 'Liberation Sans'
DPI = 300

# deck palette (deck/build_slides.js)
NAVY, BLUE, GREEN, RED = '#1F3864', '#0571B0', '#008837', '#B2182B'
SIG, OCHRE, MIDBLUE = '#FCEFD2', '#B07C1E', '#4F81BD'     # the deck palette (tmp/ppt/hla_assoc): highlight, platforms
INK, GREY, RULE, BAND, PALE = '#212121', '#595959', '#BFBFBF', '#DEEAF6', '#F4F7FB'
FLAGF, OKF, NAF = '#F6D5D8', '#E3EEF7', '#F2F2F2'

# the experiments: key, label, the registry checks each one feeds (evidence-matrix columns), the verb of its factor
EXPERIMENTS = [
    ('Q1', 'Depth: read down-sampling', ['depth'], 'reduces depth'),
    ('Q2', 'Platform: platform strata', ['platform_Q', 'platform_RE', 'lopo'], 'stratifies the cases by platform'),
    ('Q3', 'Technology: SNP arrays', ['array'], 'replaces the genotyping technology'),
    ('Q4', 'Control source: source strata', ['sources_Q', 'sources_loo'], 'stratifies the controls by source'),
]
_known = {c for _, _, cs, _ in EXPERIMENTS for c in cs} | {'discovery'}
if set(CK.RULES) != _known:
    sys.exit(f'the experiments do not cover the registry checks: {sorted(set(CK.RULES) ^ _known)}')


# ── small readers ──────────────────────────────────────────────────────────────────────────────────────────────
def tsv(*p, **kw):
    return pd.read_csv(RES.joinpath(*p), sep='\t', **kw)


def js(*p):
    return json.loads(RES.joinpath(*p).read_text())


def rtsv(name, **kw):
    """a report table (report_file)"""
    return pd.read_csv(report_file(name), sep='\t', **kw)


def rjs(name):
    return json.loads(report_file(name).read_text())


def imputation(text):
    """the imputation record of the arrays, as 11.array/array.summary.json gives it: '<tool> v<version> (<record>)
    against <panel>, ...'. The deck and docs/DESIGN.md take the tool from here only (one stats.json key)."""
    m = re.match(r'(\S+) v([\d.]+) \(([^)]*)\) against (.+?), which', text)
    if not m:
        sys.exit(f'unexpected imputation record: {text}')
    tool, version, record, panel = m.groups()
    record = re.sub(r', \d{4}-\d\d-\d\d$', '', record)          # the file date stays in the results
    return OrderedDict(tool=tool, version=version, record=record, panel=panel)


def _dur(s):
    """a trace duration ('1h 2m 3s', '250ms') in seconds"""
    if not isinstance(s, str) or s in ('-', ''):
        return 0.0
    unit = {'ms': 1e-3, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}
    return sum(float(v) * unit[u] for v, u in re.findall(r'([\d.]+)\s*(ms|d|h|m|s)\b', s))


def cores_per_process():
    """the cores each process is given (nextflow.config: clusterOptions --rsc c=; the default for the rest)"""
    cfg = (ROOT / 'nextflow.config').read_text()
    params = {k: int(v) for k, v in re.findall(r'^params\.(\w+)\s*=\s*(\d+)\s*$', cfg, re.M)}
    default = int(re.search(r"process \{.*?clusterOptions = '--rsc p=\d+:t=\d+:c=(\d+)", cfg, re.S).group(1))
    rules = []
    for pat, body in re.findall(r"withName: '([^']+)' \{(.*?)\n    \}", cfg, re.S):
        v = re.search(r'c=(\$\{[^}]+\}|\d+)', body).group(1)
        if v.startswith('${'):
            v = eval(re.sub(r'params\.(\w+)', lambda m: str(params[m.group(1)]), v[2:-1]), {'__builtins__': {}})
        rules.append((pat, int(v)))

    def cores(name):
        c = default
        for pat, v in rules:
            if re.fullmatch(pat, name):
                c = v
        return c
    return cores


def compute_used():
    """core-hours from every trace: allocated = realtime x the cores of the process; used = realtime x %cpu.
    'published' counts each task behind the published results once: the tasks the latest full run (gates,
    down-sampling and platform stages, cached or completed) and the latest report run used, each at the run that
    executed it. Tasks of the superseded R = 10 arms are therefore not counted; 'executed' counts every task of every
    run, pilots included. Undocumented processes are left out of both."""
    cores = cores_per_process()
    rows = []
    for f in sorted(RES.glob('_run_info/*/trace.txt')):
        t = pd.read_csv(f, sep='\t', dtype=str)
        t['run'] = f.parent.name
        rows.append(t)
    T = pd.concat(rows, ignore_index=True)
    T['process'] = T['name'].str.split(' ').str[0]
    T = T[~T.process.isin(UNDOCUMENTED)]
    T['core_h'] = T['realtime'].map(_dur) * T['process'].map(cores) / 3600
    T['cpu_h'] = T['realtime'].map(_dur) * pd.to_numeric(T['%cpu'].str.rstrip('%'), errors='coerce').fillna(0) / 360000
    executed = T[T.status.isin(['COMPLETED', 'FAILED'])]
    full_tag = CK.STAGE['downsampling'].tag                     # 'full': the gates, downsampling and platform stages
    report_tag = CK.STAGE['report'].tag

    def latest(tag):
        runs = sorted({r for r in T.run if r.split('.')[0] == tag and r.split('.')[-1].isdigit()},
                      key=lambda r: int(r.split('.')[1]))
        return runs[-1:]
    final_runs = latest(full_tag) + latest(report_tag)
    hashes = set(T[T.run.isin(final_runs)].hash)
    done = T[T.status == 'COMPLETED'].drop_duplicates('hash', keep='last')
    pub = done[done.hash.isin(hashes)]
    stage_of = {p.name: p.stage for p in PROCESSES}
    by_stage = pub.groupby(pub.process.map(stage_of)).core_h.sum()
    return OrderedDict(
        runs=int(T.run.nunique()), executed_tasks=int(len(executed)), executed_core_h=float(executed.core_h.sum()),
        executed_cpu_h=float(executed.cpu_h.sum()), final_runs=final_runs, tasks=int(len(pub)),
        processes=int(pub.process.nunique()), core_h=float(pub.core_h.sum()), cpu_h=float(pub.cpu_h.sum()),
        core_h_by_stage={s.key: float(by_stage.get(s.key, 0.0)) for s in CK.STAGES if s.key in set(stage_of.values())},
        core_h_by_process={k: float(v) for k, v in pub.groupby('process').core_h.sum().sort_values(ascending=False).items()},
        tasks_by_process={k: int(v) for k, v in pub.groupby('process').size().items()})


# ── the numbers ────────────────────────────────────────────────────────────────────────────────────────────────
def numbers():
    N = OrderedDict()
    N['sources_of_numbers'] = OrderedDict(results=str(RES), tables=str(TABLES) if TABLES else None)
    cnt = tsv('00.samples', 'counts.tsv')
    N['cases'] = {c: int(cnt[(cnt.cohort == c) & (cnt.status == 'case')].n.sum()) for c in COHORTS}
    N['controls'] = {c: int(cnt[(cnt.cohort == c) & (cnt.status == 'control')].n.sum()) for c in COHORTS}
    N['genomes'] = {c: N['cases'][c] + N['controls'][c] for c in COHORTS}
    N['case_platform'] = {c: {p: int(cnt[(cnt.cohort == c) & (cnt.platform == p)].n.sum()) for p in CASE_PLATFORMS}
                          for c in COHORTS}
    cc = rtsv('controls.counts.tsv')
    N['source'] = {c: {r.source: int(r.controls_random) for r in
                       cc[(cc.cohort == c) & cc.subset.str.startswith('ctrlonly')].itertuples()} for c in COHORTS}
    pdp = tsv('01.depth', 'platform_depth.tsv')
    N['depth'] = {r.platform: {'n': int(r.n), 'readlen': float(r.readlen_median), 'depth': float(r.depth_median),
                               'q25': float(r.depth_q25), 'q75': float(r.depth_q75)}
                  for r in pdp[pdp.cohort == 'full_mainland'].itertuples()}
    dsm = js('01.depth', 'depth_summary.json')
    N['target_depth'] = float(dsm['target_depth'])
    N['fraction_median'] = float(dsm['fraction_median'])
    N['fraction_by_platform'] = {k: [float(x) for x in v] for k, v in dsm['fraction_by_platform'].items()}
    ss = js('00.samples', 'sample_summary.json')
    N['deep'] = {'n': int(ss['deep_cases']), 'by_platform': {k: int(v) for k, v in ss['deep_cases_by_platform'].items()},
                 'share_of_cases': {c: ss['deep_cases'] / N['cases'][c] for c in COHORTS}}

    # the down-sampled cases (user, 2026-11-03): every case platform but DNBseq-G400RS 15x to the controls' median;
    # the deep cases each to it, the DNBSeq-T7 cases by one platform fraction (T7_FRACTIONS -> 01.depth/fractions.all)
    fr = tsv('01.depth', 'fractions.all.tsv')
    t7 = fr[fr.platform == 'DNBSeq-T7 30x']
    st = tsv('00.samples', 'sample_table.tsv', dtype=str).set_index('sample_id')
    ds_ids = set(tsv('03.calls', 'downsampling', 'deep_cases.tsv', dtype=str).sample_id)
    deep_ids = set(fr[fr.platform.isin(N['deep']['by_platform'])].sample_id)
    by_c = {}
    for c in COHORTS:
        inc = {s for s in ds_ids if st.loc[s, f'in_{c}'] == '1'}
        by_c[c] = {'all': len(inc), 'deep': len(inc & deep_ids), 't7': len(inc - deep_ids)}
    N['downsampled'] = OrderedDict(
        n=len(ds_ids), deep=len(ds_ids & deep_ids), t7=len(ds_ids - deep_ids), by_cohort=by_c,
        share_of_cases={c: by_c[c]['all'] / N['cases'][c] for c in COHORTS},
        t7_fraction=float(t7.fraction.astype(float).iloc[0]), t7_fractions=int(t7.fraction.nunique()),
        t7_median=N['depth']['DNBSeq-T7 30x']['depth'],
        untouched=N['case_platform']['full_mainland']['DNBseq-G400RS 15x'],
        platforms=sorted(fr.platform.unique().tolist()))

    # down-sampling arms (arm plan), joint call, HLA re-typing, swaps, re-test arms
    ap = tsv('03.calls', 'downsampling', 'arm_plan.tsv')
    deep, fid = ap[ap.role == 'deep'], ap[ap.role == 'fidelity']
    reps = sorted(a for a in deep.arm.unique() if re.fullmatch(r'r\d+', a))
    og = tsv('03.calls', 'downsampling', 'orig_gvcfs.tsv')
    swaps = glob.glob(str(RES / '04.genotypes' / 'downsampling' / '*' / '*' / 'swap.*.json'))
    N['arms'] = OrderedDict(
        replicates=len(reps), replicate_names=reps, arms=int(deep.arm.nunique()), deep_genome_arms=int(len(deep)),
        fidelity_controls=int(fid.sample_id.nunique()), fidelity_platforms=sorted(fid.platform.unique().tolist()),
        called_genome_arms=int(len(ap)), original_gvcfs=int(len(og)), joint_gvcfs=int(len(ap) + len(og)),
        hla_typings=int(len(deep)), swap_records=len(swaps),
        retest_arms=1 + len([d for d in (RES / '06.retest' / 'downsampling').iterdir()
                             if d.is_dir() and (d.name == 'baseline' or d.name in reps)]))
    tq = js('02.targets', 'targets_qc.json')
    N['targets'] = OrderedDict(fixed_sites=int(tq['fixed_sites']), by_kind={k: int(v) for k, v in tq['fixed_by_kind'].items()},
                               intervals=int(tq['hc_intervals']), interval_kb=tq['hc_bp'] / 1000,
                               null_snps=int(tq['null_snps']), gene_variants=int(tq['gene_tested_variants_distinct']),
                               flank_kb=tq['flank'] / 1000)

    # the signals and the markers each re-test reads
    sig = tsv('02.targets', 'signal_variants.tsv', dtype=str)
    gw = sig[sig.kind == 'gwas']
    N['signals'] = OrderedDict(
        loci=int(gw.unit.nunique()), gwas_variants=int(gw.variant_id.nunique()),
        leads=int((gw.role == 'Lead').sum()), alt_leads=int((gw.role == 'Alt lead').sum()),
        cs_members=int((gw.role == 'CS member').sum()), genes=sorted(sig[sig.kind == 'gene'].unit.unique().tolist()),
        gene_sets=int(sig[sig.kind == 'gene'][['unit', 'stratum']].drop_duplicates().shape[0]),
        hla_alleles={c: int(len(tsv('06.retest', 'original', c, 'assoc.random.allele.tsv'))) for c in COHORTS},
        hla_residues={c: int(len(tsv('06.retest', 'original', c, 'assoc.random.residue.tsv'))) for c in COHORTS},
        positions={c: int(n) for c, n in tsv('06.retest', 'original', '_omnibus', 'omnibus.residue.tsv')
                   .groupby('cohort').size().items()},
        gwas_retest_variants={c: int(len(tsv('06.retest', 'original', c, 'gwas.best_guess.tsv'))) for c in COHORTS})

    # platform and control-source subsets: null models, step-2 runs, gene runs, omnibus runs
    subs = {c: sorted(d.name for d in (RES / '08.subsets' / c / 'null').iterdir() if d.is_dir()) for c in COHORTS}
    plat = {c: [s for s in subs[c] if not s.startswith('ctrl')] for c in COHORTS}
    srcs = {c: [s for s in subs[c] if s.startswith('ctrl')] for c in COHORTS}
    step2 = {c: sorted(glob.glob(str(RES / '09.platform_assoc' / c / '*.saige.tsv'))) for c in COHORTS}
    rows = {c: [sum(1 for _ in open(f)) - 1 for f in step2[c]] for c in COHORTS}
    N['subsets'] = OrderedDict(
        platform=len(plat['full_mainland']), sources=len(srcs['full_mainland']), per_cohort=len(subs['full_mainland']),
        null_models=sum(len(v) for v in subs.values()), step2_runs=sum(len(v) for v in step2.values()),
        markers_per_run={c: [min(rows[c]), max(rows[c])] for c in COHORTS}, step2_tests=int(sum(sum(v) for v in rows.values())),
        gene_runs=sum(len(glob.glob(str(RES / '09.platform_assoc' / c / 'gene.*.tsv'))) for c in COHORTS),
        omnibus_runs=len(glob.glob(str(RES / '09.platform_assoc' / '_omnibus' / 'omnibus.*.tsv'))),
        subset_counts=tsv('08.subsets', 'full_mainland', 'subsets', 'subset_counts.tsv')
        .set_index('subset')[['cases_random', 'controls_random']].to_dict('index'),
        source_subset_counts=tsv('08.subsets', 'full_mainland', 'control_sources', 'subsets', 'subset_counts.tsv')
        .set_index('subset')[['cases_random', 'controls_random']].to_dict('index'))
    N['platform_null'] = {c: js('10.platform', c, 'platform.summary.json')['null_share_p_Q_lt_0_05'] for c in COHORTS}
    rho = pd.concat([tsv('10.platform', c, 'platform.rho.tsv').assign(cohort=c) for c in COHORTS])
    N['rho'] = {'analytic': [float(rho.rho_analytic.min()), float(rho.rho_analytic.max())],
                'null_snps': [float(rho.r_null_snps.min()), float(rho.r_null_snps.max())]}
    lam = rtsv('controls.null_lambda.tsv')
    N['sources_lambda'] = [float(lam['lambda'].min()), float(lam['lambda'].max())]
    N['sources_share_z3'] = [float(lam.share_abs_z_gt3.min()), float(lam.share_abs_z_gt3.max())]
    N['share_z3_null'] = math.erfc(3 / math.sqrt(2))          # P(|z| > 3) under the null

    # SNP arrays: check/array_crosscheck's comparison re-run in each cohort (array_cohorts.py); gate G10
    asum = rjs('array.summary.json')
    asig = rtsv('array.signals.tsv')
    keys = ('shared_sites', 'qc_sites', 'concordance_cases_T7_vs_ASA24v1', 'concordance_controls_HiSeqX_vs_JSA24v1',
            'lambda_wgs', 'lambda_array', 'r_z_wgs_array', 'slope_z_array_on_wgs', 'abs_dz_gt4', 'sites_wgs_p_lt_1e-4',
            'of_which_array_same_direction_p_lt_0.01', 'signals_by_verdict')
    N['array'] = OrderedDict(
        (c, OrderedDict(cases=int(v['people']['cases']), controls=int(v['people']['controls']),
                        people=int(v['people']['cases'] + v['people']['controls']),
                        signals=int((asig.cohort == c).sum()), **{k: v[k] for k in keys},
                        # the array SNPs of the DQB1:57:V prediction, the cohort's own compare step
                        v57_features=int(js('11.array', 'cohorts', c, '04.compare', 'summary.json')
                                         ['hla_dqb1_57V']['features'])))
        for c, v in asum['cohorts'].items())
    g10 = asum['gate_G10']
    N['G10'] = OrderedDict(passed=bool(g10['pass']), tables=len(g10['tables']),
                           max_rel_diff=max(float(t['max_rel_diff']) for t in g10['tables'].values()))
    # the imputation of the arrays: the one record of the tool (deck: stats.json key 'imputation'). To name another
    # tool, change this one line.
    N['imputation'] = imputation(asum['imputation'])

    # follow-up: the re-applied genotype filters, per cohort
    gf = rtsv('followup.gene_filters.tsv')
    N['followup'] = OrderedDict((c, OrderedDict(gene_variants=int((gf.cohort == c).sum()),
                                                genomes=int(gf[gf.cohort == c].samples.iloc[0]),
                                                agreement=float(gf[gf.cohort == c].genotype_agreement.min())))
                                for c in COHORTS)

    # the evidence matrix and the decision
    lg = rtsv('evidence.long.tsv')
    ev = rtsv('evidence.signals.tsv')
    chk = lg[lg.check != 'discovery']
    N['matrix'] = OrderedDict(
        signals=int(len(ev)), by_kind={k: int(v) for k, v in ev.kind.value_counts().items()},
        checks=int(len([c for c in ev.columns if c in CK.RULES and c != 'discovery'])), cohorts=len(COHORTS),
        readings=int(len(chk)), discovery_readings=int((lg.check == 'discovery').sum()),
        not_applicable=int((chk.reading == 'not applicable').sum()),
        applying=int((chk.reading != 'not applicable').sum()),
        flag_cells=int((chk.reading == 'flag').sum()),
        flag_signals=sorted(ev[ev['flags'] != 'none'].signal.tolist()),
        readings_by_kind={k: int(v) for k, v in chk.reading.value_counts().items()})

    # the registry (documented processes only)
    N['registry'] = OrderedDict(processes=len(PROCESSES), gates=len(CK.GATES), checks=len(CK.CHECKS),
                                stages={s: len([p for p in PROCESSES if p.stage == s])
                                        for s in ['inputs', 'gates', 'downsampling', 'platform', 'report']},
                                min_cohorts=CK.MIN_COHORTS, loo_z_expected=CK.LOO_Z_EXPECTED,
                                loo_z_observed=CK.LOO_Z_OBSERVED, depth_dz=CK.DEPTH_DZ, q_alpha=CK.Q_ALPHA,
                                q_min_maf=CK.Q_MIN_MAF)
    # the fidelity gates (results/_gates; G10 in the array summary)
    g = OrderedDict()
    g['G1'] = {c: js('_gates', f'G1.{c}.json') for c in COHORTS}
    g['G2'] = {f'{c}|{st_}': js('_gates', f'G2.{c}.{st_}.json') for c in COHORTS
               for st_ in ('moderate_high', 'low_moderate_high')}
    g3 = js('_gates', 'G3.json')
    g['G3'] = {'files': len(g3['files']), 'identical': sum(v['published'] == v['reproduced'] for v in g3['files'].values())}
    g4 = js('_gates', 'G4.json')
    g['G4'] = {'single_rows': g4['single_markers']['rows_ours'], 'omnibus_rows': g4['omnibus']['rows_ours'],
               'worst': max(v for k, v in list(g4['single_markers'].items()) + list(g4['omnibus'].items())
                            if k.startswith('max_rel_diff'))}
    g['G5'] = {c: js('_gates', f'G5.{c}.json') for c in COHORTS}
    g['G6'] = {f"{r['comparison']}|{r['role']}": r for r in js('_gates', 'downsampling', 'G6.json')['rows']}
    g['G8'] = js('_gates', 'downsampling', 'G8.json')
    g['G9'] = js('_gates', 'downsampling', 'G9.json')
    g['G10'] = N['G10']
    g['pass'] = all(js('_gates', f.name)['pass'] for f in sorted((RES / '_gates').glob('G*.json'))) and \
        all(js('_gates', 'downsampling', f'{k}.json')['pass'] for k in ('G6', 'G8', 'G9')) and N['G10']['passed']
    N['gates'] = g
    # the emulation's QC profile (positive control of Q1), per cohort: the median share of genotypes failing any QC
    # rule; the cohort's own pre-QC genotypes by group, the re-called arms per platform of the down-sampled cases
    qr = rtsv('ds.qc_reasons.tsv')
    q = qr[qr.arm.isin(['original', 'baseline'] + reps)]

    def share(c, source, group, platform, arms):
        x = q[(q.cohort == c) & (q.source == source) & (q.group == group) & (q.platform == platform) & q.arm.isin(arms)]
        return [float(v) for v in x['any']] if len(x) else []
    N['qc_fail'] = OrderedDict()
    for c in COHORTS:
        d = OrderedDict(controls=share(c, 'cohort', 'control', 'all', ['original']),
                        deep_original=share(c, 'cohort', 'deep', 'all', ['original']),
                        t7_original=share(c, 'cohort', 'DNBSeq-T7 30x', 'all', ['original']),
                        g400rs15_original=share(c, 'cohort', 'DNBseq-G400RS 15x', 'all', ['original']))
        for p in ['NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBSeq-T7 30x']:
            d[f'baseline|{p}'] = share(c, 'recalled', 'deep', p, ['baseline'])
            d[f'reps|{p}'] = share(c, 'recalled', 'deep', p, reps)
        d['reps|all'] = share(c, 'recalled', 'deep', 'all', reps)
        N['qc_fail'][c] = OrderedDict((k, (v[0] if len(v) == 1 and not k.startswith('reps') else v)) for k, v in d.items())
    # the replicate count: the user's decision for every platform (three); the pilot's measured costs beside it
    rd = js('_pilot', 'R_decision.json')
    N['pilot'] = {'R_rule': rd['R'], 'rule': rd['rule'], 'projection': rd['projection'],
                  'core_h_per_task': rd['core_h_per_task'], 'median_realtime_s': rd['median_realtime_s'],
                  'tasks_measured': rd['tasks_measured'], 'R': len(reps)}
    met = lg[(lg.check == 'platform_Q') & (lg.signal == 'FTO') & (lg.cohort == 'full_mainland')].metric.iloc[0]
    N['q_threshold'] = float(re.search(r'threshold ([\d.]+)', met).group(1))
    lo = rtsv('evidence.loo.tsv')
    lp = lo[(lo.check == 'lopo') & (lo.kind != 'gene')]
    N['loo_sd'] = {'t7': [float(lp[lp.left_out == 't7'].sd_conditional.min()), float(lp[lp.left_out == 't7'].sd_conditional.max())],
                   'small': [float(lp[lp.left_out != 't7'].sd_conditional.min()), float(lp[lp.left_out != 't7'].sd_conditional.max())]}
    N['compute'] = compute_used()
    return N


# ── the process graph ──────────────────────────────────────────────────────────────────────────────────────────
def _statements(block):
    """the workflow block as statements: a line continues while brackets are open or the next line starts with '.'"""
    lines = [re.sub(r'//.*$', '', ln).rstrip() for ln in block.splitlines()]
    # block structure (workflow {, if (...) {, } else {, a block's closing }) is not part of a statement; a lone }
    # inside an open statement closes a closure
    opener = re.compile(r'^\s*(?:workflow\s*\{|(?:\}\s*else\s*)?if\s*\(.*\)\s*\{|\}\s*else\s*\{)\s*$')
    out, buf, depth = [], '', 0
    for i, ln in enumerate(lines):
        if not ln.strip() or (depth <= 0 and (opener.match(ln) or ln.strip() == '}')):
            continue
        buf += ' ' + ln.strip()
        depth += sum(ln.count(o) for o in '([{') - sum(ln.count(c) for c in ')]}')
        nxt = next((x.strip() for x in lines[i + 1:] if x.strip()), '')
        if depth <= 0 and not nxt.startswith('.'):
            out.append(buf.strip())
            buf, depth = '', 0
    return out


def _call_args(text, start):
    """the text inside the parentheses that open at text[start]"""
    depth = 0
    for j in range(start, len(text)):
        depth += text[j] == '('
        depth -= text[j] == ')'
        if depth == 0:
            return text[start + 1:j]
    return text[start + 1:]


def channel_edges():
    """(source, target): an output channel of source is an input of target in platform_depth.nf (every stage; the
    report processes read the published folders, see folder_reads, and pass their tables on by channel)"""
    names = [p.name for p in PROCESSES]
    nf = (ROOT / 'platform_depth.nf').read_text()
    block = nf[nf.index('\nworkflow {'):]
    var_src = defaultdict(set)
    edges = set()
    word = re.compile(r'(?<![.\w])([A-Za-z_]\w*)\b')      # a name, not an attribute after a dot
    for st in _statements(block):
        m = re.match(r'(?:def\s+)?([a-z_]\w*)\s*=(?!=)\s*(.*)$', st)
        lhs, rhs = (m.group(1), m.group(2)) if m else (None, st)
        called = set()
        for cm in re.finditer(r'\b(' + '|'.join(names) + r')\s*\(', rhs):
            for w in set(word.findall(_call_args(rhs, cm.end() - 1))):
                edges.update((s, cm.group(1)) for s in var_src.get(w, ()))
            called.add(cm.group(1))
        if lhs:
            refs = set().union(*[var_src.get(w, set()) for w in set(word.findall(rhs))])
            var_src[lhs] = (called or refs) | (var_src[lhs] if lhs in word.findall(rhs) else set())
    return sorted(edges)


def _folder(path):
    """the top folder of a results path ('06.retest/<arm>/...' -> '06.retest')"""
    return path.strip().split('/')[0]


def published_folders():
    """{process: the top folders under results/ it publishes} (registry: publishes)"""
    return {p.name: list(dict.fromkeys(_folder(x) for x in p.publishes.split(',') if x.strip()))
            for p in PROCESSES}


def folder_reads():
    """(folder, report process, the path as the registry's `reads` gives it): the published folders each
    report-stage process reads"""
    out = []
    ours = {f for fs in published_folders().values() for f in fs}
    for p in PROCESSES:
        if p.stage != 'report':
            continue
        for tok in re.findall(r'\b(\d\d\.[\w./]+)', p.reads):
            if _folder(tok) in ours:                           # results/ folders, not another component's
                out.append((_folder(tok), p.name, tok.rstrip('/.')))
    return out


def registry_mentions():
    """(upstream, process) for every process the registry's `reads` names explicitly"""
    names = [p.name for p in PROCESSES]
    out = set()
    for p in PROCESSES:
        for q in names:
            if q != p.name and re.search(r'\b' + q + r'\b', p.reads):
                out.add((q, p.name))
    return out


# ── rendering ──────────────────────────────────────────────────────────────────────────────────────────────────
MANIFEST = {}


def _esc(t):
    return str(t).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _png_size(path):
    b = Path(path).read_bytes()
    return int.from_bytes(b[16:20], 'big'), int.from_bytes(b[20:24], 'big')


def render_dot(name, src, min_font, svg=False):
    """write docs/design/<name>.dot, render <name>.png at DPI (and .svg) with Liberation Sans; record the size"""
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f'{name}.dot'
    f.write_text(src)
    env = dict(os.environ, FONTCONFIG_FILE=str(FONTS_CONF))
    subprocess.run([DOT, f'-Gdpi={DPI}', '-Tpng', '-o', str(OUT / f'{name}.png'), str(f)], check=True, env=env)
    if svg:
        subprocess.run([DOT, '-Tsvg', '-o', str(OUT / f'{name}.svg'), str(f)], check=True, env=env)
    w, h = _png_size(OUT / f'{name}.png')
    MANIFEST[name] = {'png': f'docs/design/{name}.png', 'source': f'docs/design/{name}.dot', 'tool': 'graphviz dot',
                      'width_in': round(w / DPI, 3), 'height_in': round(h / DPI, 3), 'min_font_pt': min_font}


def fmt_int(x):
    return f'{int(round(x)):,}'


def rng(xs, f=fmt_int, sep='\u2013'):
    lo, hi = min(xs), max(xs)
    return f(lo) if f(lo) == f(hi) else f'{f(lo)}{sep}{f(hi)}'


def three(d, f=fmt_int):
    """N / I / F values of a per-cohort dict"""
    return ' / '.join(f(d[c]) for c in COHORTS)


# ── the process DAG ────────────────────────────────────────────────────────────────────────────────────────────
STAGE_TINT = {'inputs': '#F2F2F2', 'gates': '#E3EEF7', 'downsampling': '#FBEFF0', 'platform': '#EAF4EE',
              'report': '#F4F7FB'}
STAGE_TITLE = {'inputs': 'Inputs', 'gates': 'Gates: the original arm re-tested', 'downsampling': 'Read down-sampling arms',
               'platform': 'Platform and control-source subsets', 'report': 'Report (never cached)'}


BROADCAST = ['SAMPLE_TABLE', 'BUILD_TARGETS']          # read by most stages: drawn to the stage, not to each process


def _stage_of(name):
    return next(p.stage for p in PROCESSES if p.name == name)


def process_dag(part='full'):
    """the processes by stage, gates marked (green), channel edges (solid), the published folders the report stage
    reads (folder nodes; dashed). part 'full' (docs): every process and edge, name, scripts and published folders.
    Parts 'A' (inputs, gates, down-sampling) and 'B' (platform, report) for the deck: one line per process; the
    sample table and targets drawn to each stage they feed; upstream processes of the other part as grey stubs."""
    E = channel_edges()
    pubs = published_folders()
    reads = folder_reads()
    read_folders = sorted({f for f, _, _ in reads})
    gates = defaultdict(list)
    for g in CK.GATES:
        gates[g.process].append(g.key)
    stages = {'full': ['inputs', 'gates', 'downsampling', 'platform', 'report'],
              'A': ['inputs', 'gates', 'downsampling'], 'B': ['platform', 'report']}[part]
    shown = {p.name for p in PROCESSES if p.stage in stages}
    compact = part != 'full'
    fs = 9 if compact else 10
    L = ['digraph process_dag {',
         f'  graph [rankdir={"LR" if part == "full" else "TB"}, newrank=true, compound=true, fontname="{FONT}", '
         f'fontsize={fs + 0.5}, nodesep={0.07 if compact else 0.10}, ranksep={0.11 if compact else 0.35}, pad=0.04, '
         f'splines=true, outputorder=edgesfirst];',
         f'  node [shape=box, style="rounded,filled", fillcolor="white", color="{NAVY}", penwidth=0.9, '
         f'fontname="{FONT}", fontsize={fs}, fontcolor="{INK}", margin="0.05,0.012", height=0.16];',
         f'  edge [color="{NAVY}", penwidth=0.7, arrowsize=0.45];']
    for st in stages:
        procs = [p for p in PROCESSES if p.stage == st]
        tag = CK.STAGE['downsampling'].tag if st == 'platform' else CK.STAGE[st].tag
        L.append(f'  subgraph cluster_{st} {{ label=<<B>{_esc(STAGE_TITLE[st])}</B>  ({len(procs)} processes; run '
                 f'tag {tag}.N)>; labeljust=l; style="rounded,filled"; fillcolor="{STAGE_TINT[st]}"; color="{RULE}"; '
                 f'fontcolor="{NAVY}";')
        for p in procs:
            g = gates.get(p.name, [])
            badge = f'  <FONT COLOR="{NAVY}"><B>{" ".join(g)}</B></FONT>' if g else ''
            if compact and p.stage == 'report' and any(r[1] == p.name for r in reads):
                fl = list(dict.fromkeys(f + '/' for f, q, _ in reads if q == p.name))
                lab = (f'<{_esc(p.name)}{badge}<BR/><FONT COLOR="{GREY}">reads ' +
                       '<BR/>'.join(_esc(', '.join(fl[i:i + 2])) for i in range(0, len(fl), 2)) + '</FONT>>')
            elif compact:
                lab = f'<{_esc(p.name)}{badge}>'
            else:
                sub = ', '.join(f + '/' for f in pubs[p.name]) if p.publishes else 'channel only'
                lab = (f'<<B>{_esc(p.name)}</B>{badge}<BR/><FONT POINT-SIZE="{fs - 1}" COLOR="{GREY}">'
                       f'{_esc(p.scripts)}  \u2192 {_esc(sub)}</FONT>>')
            extra = f', color="{NAVY}", penwidth=1.8' if g else ''
            L.append(f'    "{p.name}" [label={lab}{extra}];')
        L.append('  }')
    stub = f'style="rounded,dashed", color="{GREY}", fontcolor="{GREY}"'
    if part == 'full':
        L.append(f'  subgraph cluster_results {{ label=<<B>Published folders the report stage reads</B>>; labeljust=l; '
                 f'style="rounded,dashed"; color="{RULE}"; fontcolor="{NAVY}";')
        for f in read_folders:
            src = [p for p, fl in pubs.items() if f in fl and _stage_of(p) != 'report']
            if part == 'B' and not all(x in shown for x in src):
                rows_ = [', '.join(src[i:i + 2]) for i in range(0, len(src), 2)]
                lab = f'<<B>{f}/</B><BR/>' + '<BR/>'.join(('\u2190 ' if i == 0 else '') + _esc(r)
                                                         for i, r in enumerate(rows_)) + '>'
            else:
                lab = f'<<B>{f}/</B>>'
            L.append(f'    "dir_{f}" [label={lab}, shape=folder, fillcolor="{PALE}", color="{GREY}", fontcolor="{GREY}"];')
        L.append('    { rank=same; ' + ' '.join(f'"dir_{f}";' for f in read_folders) + ' }')
        L.append('  }')
    drawn = set()
    for a, b in E:
        if a in shown and b in shown:
            if a in BROADCAST and _stage_of(b) != _stage_of(a):
                key = (a, _stage_of(b))
                if key not in drawn:
                    first = next(p.name for p in PROCESSES if p.stage == _stage_of(b))
                    L.append(f'  "{a}" -> "{first}" [lhead=cluster_{_stage_of(b)}, color="{GREY}", style=bold];')
                    drawn.add(key)
                continue
            L.append(f'  "{a}" -> "{b}";')
        elif part == 'A' and a in shown:
            if ('A', 'out') not in drawn:
                L.append(f'  "next" [label="to part 2: platform subsets and report", {stub}];')
                drawn.add(('A', 'out'))
            if (a, 'next') not in drawn:
                L.append(f'  "{a}" -> "next" [color="{GREY}", style=dashed];')
                drawn.add((a, 'next'))
        elif part == 'B' and b in shown:
            if (a, 'stub') not in drawn:
                L.append(f'  "{a}" [label="{a} (part 1)", {stub}];')
                drawn.add((a, 'stub'))
            if a in BROADCAST:
                key = (a, _stage_of(b))
                if key not in drawn:
                    first = next(p.name for p in PROCESSES if p.stage == _stage_of(b))
                    L.append(f'  "{a}" -> "{first}" [lhead=cluster_{_stage_of(b)}, color="{GREY}", style=bold];')
                    drawn.add(key)
                continue
            L.append(f'  "{a}" -> "{b}" [color="{GREY}"];')
    if part == 'B':                     # the report reads what the platform stage published: below it
        L.append(f'  "PLATFORM_STATS" -> "{next(p.name for p in PROCESSES if p.stage == "report")}" [style=invis];')
    if part == 'full':
        for p, fl in pubs.items():
            for f in fl:
                if f in read_folders and _stage_of(p) != 'report' and p in shown:
                    L.append(f'  "{p}" -> "dir_{f}" [color="{GREY}", arrowhead=none, penwidth=0.6];')
        for f, p, tok in sorted(set(reads)):
            L.append(f'  "dir_{f}" -> "{p}" [style=dashed, color="{GREY}"];')
    L.append('}')
    render_dot({'full': 'process_dag', 'A': 'process_dag_part1', 'B': 'process_dag_part2'}[part],
               '\n'.join(L) + '\n', fs, svg=part == 'full')


# ── master workflow and numbers flow (graphviz) ────────────────────────────────────────────────────────────────
GENES_I = lambda gs: ', '.join(f'<I>{g}</I>' for g in gs)          # noqa: E731  gene names in italics


def _wrap(t, width):
    """wrap plain text to lines (HTML-escaped); markup-free text only"""
    return [_esc(x) for x in textwrap.wrap(t, width)] or ['']


def box(title, lines, fs=9.5, tfs=10.5, head=NAVY, head_fc='white', body='white', w=None, border=None, bh=None):
    """an HTML-like label in the deck's card style: a coloured title bar over a white body, one outer border;
    w (inches) fixes the width of both cells so that boxes of one row line up; bh (inches) the body's height, so that
    they also end on one line (the text at the top)"""
    tw = f' WIDTH="{int(w * 72)}"' if w else ''
    th = f' HEIGHT="{int(round(bh * 72))}" VALIGN="TOP"' if bh else ''
    return (f'<<TABLE BORDER="1" COLOR="{border or head}" CELLBORDER="0" CELLSPACING="0" CELLPADDING="3">'
            f'<TR><TD BGCOLOR="{head}"{tw}><FONT COLOR="{head_fc}" POINT-SIZE="{tfs}"><B>{title}</B></FONT></TD></TR>'
            f'<TR><TD BGCOLOR="{body}"{tw}{th}><FONT POINT-SIZE="{fs}">{"<BR/>".join(lines)}</FONT></TD></TR></TABLE>>')


def plainbox(lines, fs=9.5, w=None, border=NAVY):
    """a box without a title bar (a later step of a column)"""
    tw = f' WIDTH="{int(w * 72)}"' if w else ''
    return (f'<<TABLE BORDER="1" COLOR="{border}" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">'
            f'<TR><TD BGCOLOR="white"{tw}><FONT POINT-SIZE="{fs}">{"<BR/>".join(lines)}</FONT></TD></TR></TABLE>>')


def _graph_head(name, rankdir='TB', nodesep=0.12, ranksep=0.22):
    return [f'digraph {name} {{',
            f'  graph [rankdir={rankdir}, newrank=true, fontname="{FONT}", fontsize=10, nodesep={nodesep}, '
            f'ranksep={ranksep}, pad=0.04, splines=true];',
            f'  node [shape=box, style="rounded", color="{NAVY}", penwidth=1.0, fontname="{FONT}", fontsize=10, '
            f'fontcolor="{INK}", margin="0.02,0.02"];',
            f'  edge [color="{NAVY}", penwidth=1.0, arrowsize=0.6, fontname="{FONT}", fontsize=9.5, fontcolor="{GREY}"];']


def master_workflow(N):
    """data -> fidelity gates -> four experiments (what each changes, its units) -> the checks each feeds ->
    the evidence matrix -> the decision; the key numbers on it"""
    S, A, M, R = N['signals'], N['arms'], N['matrix'], N['registry']
    sub, DS, AR = N['subsets'], N['downsampled'], N['array']
    # one idea per line, no unit counts (user, 2026-11-04: the overview must not pile text; the counts are in the notes,
    # the numbers flow and the slides of each question)
    ex = {
        'Q1': ['lowers the read depth', f'{DS["n"]} cases to the controls\' {N["target_depth"]:.1f}x',
               f'baseline and {A["replicates"]} replicates'],
        'Q2': ['stratifies the cases by platform', 'each platform alone and left out', 'all controls in every subset'],
        'Q3': ['replaces the technology', 'SNP arrays on the same people', 'DNBSeq-T7 cases, HiSeq X controls'],
        'Q4': ['stratifies the controls by source', 'each source alone and left out', 'all cases in every subset'],
    }
    W = 2.4
    L = _graph_head('master_workflow', ranksep=0.17, nodesep=0.1)
    L.append('  data [label=' + box('Main study: three nested cohorts', [
        f'cases {three(N["cases"])} on four platforms',
        f'controls {three(N["controls"])} on HiSeq X'], head=GREY, w=3.6) + ', shape=plain];')
    L.append('  sig [label=' + box('Signals of the main analysis', [
        f'{S["loci"]} GWAS leads; HLA residue DQB1:57:V',
        f'Tier 1–2 genes {GENES_I(["STBD1", "OSBP", "RPH3A", "JAK2"])}'], head=GREY, w=3.6)
        + ', shape=plain];')
    L.append('  gates [label=' + box('Fidelity checks: both pass before any comparison', [
        'same data, same analysis: every main result reproduced',
        'reads re-processed, none removed: genotypes and HLA types unchanged'], head=NAVY, w=5.6)
        + ', shape=plain];')
    L.append('  { rank=same; data; sig; }')
    L.append('  data -> gates; sig -> gates;')
    # review round 8 (C-06, D-07): the readouts each experiment feeds, named as the evidence matrix names them
    READOUTS = {'Q1': ['Read down-sampling'],
                'Q2': ['Platform estimates (Q, RE)', 'Without one small platform', 'Without DNBSeq-T7 cases'],
                'Q3': ['SNP arrays'],
                'Q4': ['Source estimates (Q)', 'Without one source']}
    assert sum(len(v) for v in READOUTS.values()) == M['checks'], 'master_workflow: readouts and checks differ in number'
    # review round 13 (F-04): one body height for the four readout boxes (the longest list), so that their tops, their
    # lists and the arrows into the matrix start on one line
    RH = max(len(v) for v in READOUTS.values()) * 9.5 * 1.2 / 72 + 6 / 72
    for k, name, checks, verb in EXPERIMENTS:
        L.append(f'  {k} [label=' + box(f'{k} · {name}', [_esc(x) for x in ex[k]], w=W) + f', shape=plain, group={k}];')
        L.append(f'  r{k} [label=' + box('Readouts in the matrix', [_esc(x) for x in READOUTS[k]], head=BAND,
                                         head_fc=NAVY, body=PALE, border=RULE, w=W, tfs=9.5, bh=RH)
                 + f', shape=plain, group={k}];')
        L.append(f'  gates -> {k}; {k} -> r{k}; r{k} -> matrix;')
    L.append('  { rank=same; ' + ' '.join(k for k, *_ in EXPERIMENTS) + '; }')
    L.append('  { rank=same; ' + ' '.join('r' + k for k, *_ in EXPERIMENTS) + '; }')
    L.append('  matrix [label=' + box('Evidence matrix: every signal through every readout, in every cohort', [
        _esc(f'{sum(len(v) for v in READOUTS.values())} readouts × {M["cohorts"]} cohorts per signal, each read by a fixed rule')],
        head=NAVY, w=6.8) + ', shape=plain];')
    L.append('  decision [label=' + box('Decision: flags, never removal', [
        _esc('a flag needs its rule met in any one cohort; the group decides what to do with it')], head=GREY, w=6.8)
        + ', shape=plain];')
    L.append('  matrix -> decision;')
    L.append('}')
    render_dot('master_workflow', '\n'.join(L) + '\n', 9.5)


def numbers_flow(N):
    """CONSORT-style flow of the units: cohorts -> re-processed reads / re-cut samples / evidence re-run or counted
    -> readings -> flags; and the compute used"""
    A, T, S, SU, M, C = N['arms'], N['targets'], N['signals'], N['subsets'], N['matrix'], N['compute']
    DS, AR, FU = N['downsampled'], N['array'], N['followup']
    dpl = N['deep']['by_platform']
    hla = {c: S['hla_alleles'][c] + S['hla_residues'][c] for c in COHORTS}
    mk = [x for v in SU['markers_per_run'].values() for x in v]
    w = 3.85
    L = _graph_head('numbers_flow', ranksep=0.11, nodesep=0.18)
    L.append('  coh [label=' + box('Three nested cohorts: narrow / intermediate / full mainland', [
        _esc(f'genomes {three(N["genomes"])}; cases {three(N["cases"])}; controls {three(N["controls"])}')],
        head=GREY, w=7.2) + ', shape=plain];')
    cols = {
        'a': ('Q1 · reads re-processed', [
            [f'{DS["n"]} down-sampled cases: {DS["deep"]} deep cases (NovaSeq {dpl["NovaSeq 30x"]},',
             f'G400RS 30x {dpl["DNBSeq-G400RS 30x"]}) and {DS["t7"]} DNBSeq-T7 cases;',
             f'per cohort {three({c: DS["by_cohort"][c]["all"] for c in COHORTS})};',
             f'{A["fidelity_controls"]} HiSeq X controls as fidelity controls (f = 1)'],
            [f'{A["arms"]} arms: baseline (f = 1) and {A["replicates"]} replicates;',
             f'{fmt_int(A["deep_genome_arms"])} case genome-arms, {fmt_int(A["called_genome_arms"])} re-called in all'],
            [f'joint call: {fmt_int(A["joint_gvcfs"])} GVCFs ({A["original_gvcfs"]} of them original)',
             f'at {fmt_int(T["fixed_sites"])} fixed sites, {T["intervals"]} intervals ({T["interval_kb"]:.0f} kb)'],
            [f'{fmt_int(A["hla_typings"])} HLA re-typings; {A["swap_records"]} swap records,',
             'every other genome byte-identical'],
            [f'{A["retest_arms"]} arms × 3 cohorts re-tested: {rng(S["gwas_retest_variants"].values())} variants,',
             f'{S["gene_sets"]} gene sets, {rng(hla.values())} HLA markers, {rng(S["positions"].values())} positions']]),
        'b': ('Q2, Q4 · samples re-cut', [
            [f'{SU["per_cohort"]} subsets per cohort: {SU["platform"]} platform', f'(pooled, 4 strata, 4 left out) and '
             f'{SU["sources"]} control-', 'source (5 alone, 5 left out)'],
            [f'{SU["null_models"]} null models (3 cohorts × {SU["per_cohort"]} subsets)'],
            [f'{SU["step2_runs"]} step-2 runs of {rng(mk)} markers', f'(signal regions, {fmt_int(T["null_snps"])} null '
             'SNPs and', f'indels, HLA markers): {fmt_int(SU["step2_tests"])} tests'],
            [f'{SU["gene_runs"]} gene runs ({SU["per_cohort"]} subsets × 3 cohorts × 2 sets);',
             f'{SU["omnibus_runs"]} residue-position omnibus runs']]),
        'c': ('Q3 · evidence re-run', [
            ['Q3: people with arrays, N / I / F:',
             f'{three({c: AR[c]["people"] for c in COHORTS})} (cases {three({c: AR[c]["cases"] for c in COHORTS})});',
             f'{fmt_int(AR["full_mainland"]["shared_sites"])} shared sites; '
             f'{AR["full_mainland"]["signals"]} signals per cohort']]),
    }
    for c, (title, boxes) in cols.items():
        prev = 'coh'
        for i, lines in enumerate(boxes):
            nid = f'{c}{i}'
            lab = box(title, [_esc(x) for x in lines], w=w) if i == 0 else plainbox([_esc(x) for x in lines], w=w)
            L.append(f'  {nid} [label={lab}, shape=plain, group={c}];')
            L.append(f'  {prev} -> {nid};')
            prev = nid
        L.append(f'  {prev} -> mat{" [weight=0]" if c == "c" else ""};')
    L.append('  { rank=same; a0; b0; c0; }  a0 -> b0 -> c0 [style=invis];')
    kinds = M['by_kind']
    L.append('  mat [label=' + box('Evidence matrix and flags', [
        _esc(f'{M["signals"]} signals ({kinds["gwas"]} GWAS leads, DQB1:57:V, position DQB1 57, {kinds["gene"]} gene sets)'),
        _esc(f'× {M["checks"]} checks × 3 cohorts = {fmt_int(M["signals"] * M["checks"] * 3)} cells: '
             f'{fmt_int(M["readings"])} defined readings, {fmt_int(M["applying"])} that apply;'),
        f'<FONT COLOR="{NAVY}"><B>' + _esc(f'{M["flag_cells"]} flag readings in {len(M["flag_signals"])} signals') +
        '</B></FONT>; none removed: a flag never removes a signal'], head=NAVY, w=5.9)
        + ', shape=plain];')
    top = list(C['core_h_by_process'].items())[0] if C['core_h_by_process'] else ('-', 0.0)
    L.append('  cmp [label=' + box('Compute behind these results', [
        _esc(f'{fmt_int(C["tasks"])} tasks of {C["processes"]} processes'),
        _esc(f'{C["core_h"]:,.0f} core-hours allocated, {C["cpu_h"]:,.0f} used'),
        _esc(f'HLA re-typing ({top[0]}): {top[1]:,.0f} core-hours'),
        _esc(f'all {C["runs"]} runs with the pilots: {C["executed_core_h"]:,.0f} core-hours')], head=GREY, w=w)
        + ', shape=plain, group=c];')
    # the follow-up counts stand beside the matrix (read with each flag), fed by no question and feeding no check;
    # the connector's words are a node of their own, so the flat edge adds no rank and no width beyond them
    L.append('  fu [label=' + box('Follow-up counts', [
        _esc(f'{three({c: FU[c]["gene_variants"] for c in COHORTS})} gene variants re-filtered;'),
        _esc(f'{S["loci"]} leads by case group and source;'), 'descriptive, not a check'], head=GREY, w=3.0)
        + ', shape=plain];')
    L.append(f'  fulab [shape=none, margin="0.1,0.02", label="read with each flag", fontsize=9.5, fontcolor="{GREY}"];')
    L.append('  { rank=same; mat; fulab; fu; }')
    L.append(f'  mat -> fulab -> fu [style=dashed, dir=none, color="{GREY}"];')
    L.append('  c0 -> cmp [style=invis, minlen=3];')
    L.append('}')
    L[1] = L[1].replace('newrank=true,', 'newrank=true, compound=true,')
    render_dot('numbers_flow', '\n'.join(L) + '\n', 9.5)


# ── threat model and coverage (matplotlib, the same palette and font) ──────────────────────────────────────────
LIB_FONTS = Path('/home/b/b37974/opt/lo/root/opt/libreoffice25.8/share/fonts/truetype')


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in LIB_FONTS.glob('LiberationSans-*.ttf'):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({'font.family': FONT, 'font.size': 10, 'svg.fonttype': 'none'})
    return plt


def _refs():
    return {k: v['short'] for k, v in json.loads((ROOT / 'deck' / 'references.json').read_text()).items()}


def threats(N):
    """the threat model, one dict per way an artifact could arise: key, short label (the rows of the deck's
    blind-spot matrix), the difference, its evidence, the path to a false association, the testable prediction, and
    the experiments that test it directly"""
    R = _refs()
    d = N['depth']
    cite = lambda *ks: '; '.join(R[k] for k in ks)  # noqa: E731
    deep = [d[p]['depth'] for p in ('NovaSeq 30x', 'DNBSeq-G400RS 30x')]
    cp = N['case_platform']['full_mainland']
    n_t7, n_all = cp['DNBSeq-T7 30x'], sum(cp.values())
    # T2 any one case platform, which the strata test; T3 DNBSeq-T7 against HiSeq X: most of the cases, whose leave-out
    # the strata can test only for strong signals, so the arrays test it directly
    # (user, 2026-11-04: an error shared by every case platform is not a realistic path, NovaSeq shares the controls'
    # chemistry; genotype filters are not a threat of their own)
    rows = [
        ('T0', 'Re-test discrepancy', 'The re-test machinery itself',
         'a re-run that differs from the main analysis',
         'a difference made by re-processing, read as an artifact (or hiding one)',
         'If the machinery differed, re-running it on the original data would not reproduce the main-analysis results',
         ['G']),
        ('T1', 'Depth-dependent genotyping error',
         f'Depth: deep cases {rng(deep, lambda x: f"{x:.1f}x")}, DNBSeq-T7 cases {d["DNBSeq-T7 30x"]["depth"]:.1f}x, '
         f'controls {d[CONTROL_PLATFORM]["depth"]:.1f}x',
         cite('garner', 'kishikawa'), 'depth-dependent genotype errors and missing genotypes, mostly rare variants and indels',
         'If depth made the signal, bringing every case platform to the controls’ depth would weaken it', ['Q1']),
        ('T2', 'Platform-specific genotyping error',
         'Calls of one case platform: chemistry, read length', cite('tom', 'kim'),
         'different calls or missing genotypes in that platform’s cases',
         'If one platform’s calls made it, the signal would sit in that platform’s cases and vanish without them',
         ['Q2']),
        ('T3', 'DNBSeq-T7 genotyping error relative to HiSeq X',
         f'DNBSeq-T7: {n_t7} of {n_all} cases (F), {d["DNBSeq-T7 30x"]["readlen"]:.0f}-bp reads against '
         f'{d[CONTROL_PLATFORM]["readlen"]:.0f} bp on HiSeq X', cite('clayton', 'tom'),
         'DNBSeq-T7 calls differently from HiSeq X; with most cases on it, leaving it out keeps enough cases for strong '
         'signals only',
         'If DNBSeq-T7 calls made it, the same people on SNP arrays would not show the association', ['Q3']),
        ('T4', 'Read-based HLA typing error', 'HLA alleles typed from the reads', cite('thuesen', 'kawaguchi'),
         'a depth- or platform-dependent typing error: a false allele or residue',
         'If typing error made the HLA signal, re-typing at the controls\' depth would move it and an array-based '
         'prediction would disagree', ['Q1', 'Q3']),
        ('T5', 'Control-source heterogeneity', 'Five control sources', cite('lee', 'guo'),
         'a batch or population difference of one source carries a signal',
         'If one control source made it, that source would stand out and leaving it out would break the signal',
         ['Q4']),
    ]
    keys = ('key', 'short', 'difference', 'evidence', 'path', 'prediction', 'experiments')
    # the deck's threat table: one short phrase per cell (the long texts stay in DESIGN.md and the figure)
    brief = {
        'T0': ('Differences in software, settings, or inputs between the re-test and the main analysis alter the '
               'test statistics',
               'The unchanged re-test of the original data fails to reproduce the main-analysis results'),
        'T1': ('Genotype calls and missing genotypes depend on read depth, which differs between case groups and controls; '
               'rare variants and indels are most affected',
               'The association attenuates when the cases are down-sampled to the controls\' depth'),
        'T2': ('The chemistry or read length of one case platform alters genotype calls or missing genotypes in its cases',
               'The effect is confined to that platform\'s cases and is lost when they are excluded'),
        'T3': ('DNBSeq-T7 differs from HiSeq X in chemistry and read length and carries most cases, so excluding it '
               'leaves enough cases for strong signals only',
               'The association weakens or reverses when the same people are genotyped on SNP arrays'),
        'T4': ('Read-based HLA typing depends on depth and chemistry; a typing error misassigns an allele or residue',
               'The HLA association moves at the controls\' depth, or an array-based typing does not show it'),
        'T5': ('A batch effect or population difference in one of the five control sources shifts control allele '
               'frequencies',
               'That source\'s estimate departs from the others, and the association disappears when it is excluded')}
    out = [dict(zip(keys, r)) for r in rows]
    for d in out:
        d['path_short'], d['pred_short'] = brief[d['key']]
    return out


# which experiment tests which threat: 'tests', 'partly' (only part of the mechanism); blind otherwise; which samples
# each experiment involves is stated separately (coverage_rows). Kept here
# once for the figure, the deck and DESIGN.md.
COVERAGE = {
    'T0': {'G': 'tests'},
    'T1': {'Q1': 'tests', 'Q2': 'partly', 'Q3': 'partly'},
    'T2': {'Q1': 'partly', 'Q2': 'tests', 'Q3': 'partly'},
    'T3': {'Q2': 'partly', 'Q3': 'tests'},
    'T4': {'Q1': 'tests', 'Q2': 'partly', 'Q3': 'tests'},
    'T5': {'Q4': 'tests'},
}
EXP_NAME = {'G': f'Fidelity checks (gates G1–G{len(CK.GATES)})', **{k: n for k, n, _, _ in EXPERIMENTS}}


def _wrap_px(fig, text, width_in, size, weight='normal', style='normal'):
    """wrap text into lines no wider than width_in at this font (measured, not estimated)"""
    r = fig.canvas.get_renderer()
    words, lines, cur = text.split(), [], ''
    t = fig.text(0, 0, '', fontsize=size, fontweight=weight, fontstyle=style)
    for w_ in words:
        cand = (cur + ' ' + w_).strip()
        t.set_text(cand)
        if t.get_window_extent(renderer=r).width / fig.dpi > width_in and cur:
            lines.append(cur)
            cur = w_
        else:
            cur = cand
    t.remove()
    return lines + [cur]


def threat_model(N):
    plt = _mpl()
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    W, H = 12.6, 0.3 + 0.055 + len(threats(N)) * (0.585 + 0.055)
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis('off')
    cols = [(0.0, 3.45, 'Difference between sequencing groups'), (3.72, 2.9, 'Path to a false association'),
            (6.89, 3.6, 'Testable prediction'), (10.76, 1.84, 'Experiment built for it')]
    head_h, row_h, gap = 0.3, 0.585, 0.055
    for x, w, t in cols:
        ax.text(x + 0.02, head_h / 2, t, fontsize=11, fontweight='bold', color=NAVY, va='center', ha='left')
    rows = threats(N)

    def cell(x, y, w, h, fc, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.06', fc=fc, ec=ec, lw=0.8))

    for i, t in enumerate(rows):
        diff, ev, path, pred, exps = t['difference'], t['evidence'], t['path'], t['prediction'], t['experiments']
        y = head_h + gap + i * (row_h + gap)
        yc = y + row_h / 2
        # A: the difference (bold navy) and its evidence (grey)
        cell(cols[0][0], y, cols[0][1], row_h, BAND, BAND)
        la = _wrap_px(fig, f'{t["key"]} · {diff}', cols[0][1] - 0.16, 10, 'bold')
        lb = _wrap_px(fig, ev, cols[0][1] - 0.16, 9.5)
        n = len(la) + len(lb)
        y0 = yc - (n - 1) * 0.085
        for j, ln in enumerate(la):
            ax.text(cols[0][0] + 0.08, y0 + j * 0.17, ln, fontsize=10, fontweight='bold', color=NAVY, va='center')
        for j, ln in enumerate(lb):
            ax.text(cols[0][0] + 0.08, y0 + (len(la) + j) * 0.17, ln, fontsize=9.5, color=GREY, va='center',
                    fontstyle='italic' if ev and not ev[0].isupper() else 'normal')
        # B and C
        for (x, w, _), txt, fc in ((cols[1], path, 'white'), (cols[2], pred, 'white')):
            cell(x, y, w, row_h, fc, RULE)
            ls = _wrap_px(fig, txt, w - 0.16, 10)
            for j, ln in enumerate(ls):
                ax.text(x + 0.08, yc + (j - (len(ls) - 1) / 2) * 0.17, ln, fontsize=10, color=INK, va='center')
        # D: the experiment chips
        x, w, _ = cols[3]
        if not exps:                      # tested directly by no experiment: the partial views, outlined
            part = [k for k, v in COVERAGE[t['key']].items() if v == 'partly']
            cell(x, y, w, row_h, 'white', NAVY)
            ax.text(x + w / 2, yc, 'partly ' + ', '.join(part), fontsize=10, fontweight='bold', color=NAVY,
                    ha='center', va='center')
        ch = (row_h - 0.06 * (max(len(exps), 1) - 1)) / max(len(exps), 1)
        for j, e in enumerate(exps):
            cy = y + j * (ch + 0.06)
            fc = NAVY
            cell(x, cy, w, ch, fc, fc)
            lab = EXP_NAME[e] if e == 'G' else f'{e} \u00b7 {EXP_NAME[e]}'
            ax.text(x + w / 2, cy + ch / 2, lab, fontsize=10 if len(exps) == 1 else 9.5, fontweight='bold',
                    color='white', ha='center', va='center')
        # arrows between the columns
        for (x0, w0, _), (x1, _, _) in zip(cols[:-1], cols[1:]):
            ax.add_patch(FancyArrowPatch((x0 + w0 + 0.03, yc), (x1 - 0.03, yc), arrowstyle='-|>', mutation_scale=9,
                                         color=GREY, lw=1.0))
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / 'threat_model.png', dpi=DPI)
    plt.close(fig)
    w_, h_ = _png_size(OUT / 'threat_model.png')
    MANIFEST['threat_model'] = {'png': 'docs/design/threat_model.png', 'source': 'scripts/design_figures.py',
                                'tool': 'matplotlib', 'width_in': round(w_ / DPI, 3), 'height_in': round(h_ / DPI, 3),
                                'min_font_pt': 9.5}


# the role of every genome group in each experiment (coverage figure; the deck draws the same)
ROLE_STYLE = OrderedDict([
    ('re-processed', (MIDBLUE, 'white', 'reads re-processed at control depth')),
    ('stratified', (NAVY, 'white', 'own stratum and left out')),
    ('second technology', (OCHRE, 'white', 'genotyped by a second technology')),
    ('original', (BAND, NAVY, 're-tested or compared on the original genotypes')),
    ('shared', ('#D9D9D9', INK, 'in every subset, unchanged')),
    ('not covered', ('white', GREY, 'not covered')),
])


def coverage_rows(N):
    """(row label, the samples it changes, {group: (role, share covered)}, annotation over the case part, over the
    control part); groups = the four case platforms and the five control sources of full mainland"""
    F = 'full_mainland'
    cp = N['case_platform'][F]
    arr, DS = N['array'][F], N['downsampled']
    every = {g: ('original', 1) for g in CASE_PLATFORMS + SOURCES}
    e1 = {**{g: ('not covered', 1) for g in CASE_PLATFORMS}, **{s_: ('shared', 1) for s_ in SOURCES}}
    e1.update({p: ('re-processed', 1) for p in DS['platforms']})
    e2 = {**{g: ('stratified', 1) for g in CASE_PLATFORMS}, **{s_: ('shared', 1) for s_ in SOURCES}}
    e3 = {**{g: ('not covered', 1) for g in CASE_PLATFORMS}, **{s_: ('second technology', arr['controls'] / N['controls'][F])
                                                                for s_ in SOURCES}}
    e3['DNBSeq-T7 30x'] = ('second technology', arr['cases'] / cp['DNBSeq-T7 30x'])
    e4 = {**{g: ('shared', 1) for g in CASE_PLATFORMS}, **{s_: ('stratified', 1) for s_ in SOURCES}}
    return [
        ('Fidelity checks', 'every genome, original genotypes', every,
         'every genome re-tested on its original genotypes', ''),
        ('Q1 · Read down-sampling', f'{DS["n"]} down-sampled cases', e1,
         f'{DS["n"]} of {N["cases"][F]} brought to {N["target_depth"]:.1f}x; G400RS 15x below it', 'unchanged, the depth target'),
        ('Q2 · Platform strata', f'{N["cases"][F]} cases by platform', e2, 'each platform alone and left out',
         'every control in every subset'),
        ('Q3 · SNP arrays', f'{fmt_int(arr["people"])} people on arrays', e3,
         f'{arr["cases"]} of {cp["DNBSeq-T7 30x"]} DNBSeq-T7 cases typed', f'{fmt_int(arr["controls"])} of '
         f'{fmt_int(N["controls"][F])} typed'),
        ('Q4 · Control sources', f'{fmt_int(N["controls"][F])} controls by source', e4, 'every case in every subset',
         'each source alone and left out'),
    ]


def coverage(N):
    plt = _mpl()
    from matplotlib.patches import Rectangle
    F = 'full_mainland'
    groups = [(p, N['case_platform'][F][p], PLAT_SHORT[p]) for p in CASE_PLATFORMS] + \
             [(s_, N['source'][F][s_], SRC_LABEL[s_]) for s_ in SOURCES]
    W, x0, x1, gap_cc = 12.6, 2.55, 12.5, 0.12
    total = sum(n for _, n, _ in groups)
    k = (x1 - x0 - gap_cc) / total                      # inches per genome: one scale for every group
    xs, x = {}, x0
    for i, (g, n, _) in enumerate(groups):
        if i == len(CASE_PLATFORMS):
            x += gap_cc
        xs[g] = (x, n * k)
        x += n * k
    rows = coverage_rows(N)
    top, rh, rg = 1.38, 0.34, 0.33
    H = top + len(rows) * (rh + rg) + 0.62
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis('off')
    # group labels (staggered for the narrow case groups) and the two brackets
    cx0, cx1 = xs[CASE_PLATFORMS[0]][0], xs[CASE_PLATFORMS[-1]][0] + xs[CASE_PLATFORMS[-1]][1]
    sx0, sx1 = xs[SOURCES[0]][0], xs[SOURCES[-1]][0] + xs[SOURCES[-1]][1]
    for a, b, t in ((cx0, cx1, f'Cases, {N["cases"][F]}, by platform'), (sx0, sx1, f'Controls, {fmt_int(N["controls"][F])}, by source')):
        ax.plot([a, a, b, b], [0.36, 0.28, 0.28, 0.36], color=GREY, lw=0.8)
        ax.text((a + b) / 2, 0.22, t, ha='center', va='bottom', fontsize=10, fontweight='bold', color=NAVY)
    for i, (g, n, lab) in enumerate(groups):
        gx, gw = xs[g]
        if g in CASE_PLATFORMS:
            ty = 0.52 + 0.19 * i
            ax.plot([gx + gw / 2, gx + gw / 2], [ty + 0.07, top - 0.03], color=RULE, lw=0.6)
            ax.text(cx0 + 0.02 if i == 0 else gx + gw / 2 - 0.05, ty, f'{lab} {fmt_int(n)}', fontsize=9.5, color=INK,
                    ha='left' if i == 0 else 'left', va='center', bbox=dict(fc='white', ec='none', pad=0.4))
        else:
            ax.text(gx + gw / 2, 0.68, lab, fontsize=9.5, ha='center', va='center', color=INK)
            ax.text(gx + gw / 2, 0.88, fmt_int(n), fontsize=9.5, ha='center', va='center', color=GREY)
    for r, (label, what, roles, note_c, note_k) in enumerate(rows):
        y = top + r * (rh + rg)
        ax.text(0.05, y + rh / 2 - 0.08, label, fontsize=10, fontweight='bold', color=NAVY, va='center')
        ax.text(0.05, y + rh / 2 + 0.11, what, fontsize=9.5, color=GREY, va='center')
        for g, (role, share) in roles.items():
            gx, gw = xs[g]
            fc, tc, _ = ROLE_STYLE[role]
            ax.add_patch(Rectangle((gx, y), gw, rh, fc='white', ec=RULE, lw=0.5))
            if role == 'not covered':
                continue
            if g in SOURCES and share < 1:               # array controls: the share typed, not split by source
                ax.add_patch(Rectangle((gx, y + rh * (1 - share)), gw, rh * share, fc=fc, ec='none'))
            else:
                ax.add_patch(Rectangle((gx, y), gw * share, rh, fc=fc, ec='none'))
            ax.add_patch(Rectangle((gx, y), gw, rh, fc='none', ec='white', lw=0.6))
        note = (f'Cases: {note_c}.   Controls: {note_k}' if note_k else note_c[0].upper() + note_c[1:])
        ax.text(cx0, y + rh + 0.03, note, fontsize=9.5, color=INK, va='top', ha='left')
    # legend, two rows
    items = list(ROLE_STYLE.items())
    for j, (role, (fc, tc, txt)) in enumerate(items):
        lx = 0.05 + (j % 4) * 3.15
        ly = H - 0.42 + (j // 4) * 0.24
        ax.add_patch(Rectangle((lx, ly - 0.08), 0.22, 0.16, fc=fc, ec=RULE, lw=0.5))
        ax.text(lx + 0.28, ly, txt, fontsize=9.5, va='center', color=INK)
    fig.savefig(OUT / 'coverage.png', dpi=DPI)
    plt.close(fig)
    w_, h_ = _png_size(OUT / 'coverage.png')
    MANIFEST['coverage'] = {'png': 'docs/design/coverage.png', 'source': 'scripts/design_figures.py', 'tool': 'matplotlib',
                            'width_in': round(w_ / DPI, 3), 'height_in': round(h_ / DPI, 3), 'min_font_pt': 9.5}


# ── docs/DESIGN.md (written from here; --check-doc compares) ──────────────────────────────────────────────────────
def _pct(x, d=1):
    return f'{100 * x:.{d}f} %'


def _f(x, d=3):
    return f'{x:.{d}f}'


def _sci(x, sig=2):
    """m × 10⁻ᵉ with Unicode superscripts"""
    e = math.floor(math.log10(abs(x)))
    m = x / 10 ** e
    if round(m, sig - 1) >= 10:
        e, m = e + 1, m / 10
    sup = str(e).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return f'{m:.{sig - 1}f} × 10{sup}'


def _md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    out += ['| ' + ' | '.join(str(c).replace('|', '\\|').replace('\n', ' ') for c in r) + ' |' for r in rows]
    return out


def experiment_cards(N):
    """the four experiment cards of DESIGN.md: one fixed template (the deck's cards say the same, shorter)"""
    F, S, A, T, SU = 'full_mainland', N['signals'], N['arms'], N['targets'], N['subsets']
    d, cp, src = N['depth'], N['case_platform'], N['source']
    R, G, DS, AR, IM = N['registry'], N['gates'], N['downsampled'], N['array'], N['imputation']
    fr = N['fraction_by_platform']
    q = N['qc_fail']
    t3 = lambda f: ' / '.join(f(c) for c in COHORTS)  # noqa: E731
    cases3 = t3(lambda c: fmt_int(N['cases'][c]))
    ctrls3 = t3(lambda c: fmt_int(N['controls'][c]))
    reps = ', '.join(A['replicate_names'])

    def pr(x):
        return _pct(x, 2) if isinstance(x, float) else ('–'.join(_pct(v, 1) for v in (min(x), max(x))) if x else 'pending')
    cards = OrderedDict()
    cards['Q1'] = OrderedDict([
        ('Question', 'Does depth make the signal? If it did, bringing every case platform to the controls\' depth would '
                     'weaken it (threat T1; also T4 for the HLA typing).'),
        ('Factor varied', f'Reduced: the read depth of {DS["n"]} cases, so that every case platform but DNBseq-G400RS 15x '
                          f'has the HiSeq X controls\' median effective depth, D* = {N["target_depth"]:.2f}x. Each of the '
                          f'{DS["deep"]} deep cases (NovaSeq, DNBSeq-G400RS 30x) keeps f = D* / D_i of its read pairs '
                          f'(f {fr["NovaSeq 30x"][0]:.2f}–{fr["NovaSeq 30x"][2]:.2f} NovaSeq, '
                          f'{fr["DNBSeq-G400RS 30x"][0]:.2f}–{fr["DNBSeq-G400RS 30x"][2]:.2f} G400RS 30x); the '
                          f'{DS["t7"]} DNBSeq-T7 cases keep one platform fraction, f = D* / {DS["t7_median"]:.2f}x = '
                          f'{DS["t7_fraction"]}, which keeps their own spread of depth. Mates stay together (a hash of the '
                          f'read name and the seed). Arms: a baseline at f = 1 and {A["replicates"]} replicates ({reps}).'),
        ('Held constant', 'Every other genome (byte-identical; G7) and the variant set (the fixed sites); the cohort\'s calling, '
                          'joint genotyping and genotype QC; the null models of the main analysis, each analysis\' model and threshold; '
                          'for HLA the same HLA-HD install, dictionary and options and the original tested columns.'),
        ('Perturbed (N / I / F)', f'Down-sampled cases {t3(lambda c: fmt_int(DS["by_cohort"][c]["all"]))} '
                                  f'({t3(lambda c: _pct(DS["share_of_cases"][c], 0))} of the cases {cases3}): the '
                                  f'{DS["deep"]} deep cases in every cohort and the DNBSeq-T7 cases each cohort holds '
                                  f'({t3(lambda c: fmt_int(DS["by_cohort"][c]["t7"]))}); the {DS["untouched"]} G400RS 15x '
                                  f'cases ({d["DNBseq-G400RS 15x"]["depth"]:.1f}x) are not touched; {A["fidelity_controls"]} '
                                  'HiSeq X controls are re-processed at f = 1 only (fidelity).'),
        ('Comparator', 'Depth effect: each replicate against the baseline (identical processing, only the reads differ). '
                       'Re-processing: the baseline against the original genotypes. Total: replicate against original.'),
        ('Model and statistic', 'Each signal\'s own analysis: logistic mixed model step 2 on the null model of the main analysis (GWAS '
                                'leads, HLA alleles and residues), CMC and SKAT-O (genes), the residue-position omnibus; z from '
                                'the two-sided P.'),
        ('Readout', 'Replicate-median |z| against the original |z|; replicates still significant; genotype concordance '
                    '(confusion counts, non-reference concordance, heterozygotes retained, call rate) per platform; the share '
                    'of genotypes each QC rule no-calls, per cohort; HLA calls and dosages per arm.'),
        ('Decision rule', CK.RULES['depth'] + '; in any one cohort.'),
        ('Positive control', f'The manipulation must work: achieved / target depth over the calling intervals (G8, '
                             f'{fmt_int(G["G8"]["genome_arms"])} genome-arms, median {_f(G["G8"]["ratio_median"])}, '
                             f'{_f(G["G8"]["ratio_min"])}–{_f(G["G8"]["ratio_max"])}); the share of genotypes failing a QC rule '
                             f'moves to the controls\' (full mainland: deep cases {pr(q[F]["deep_original"])} original, '
                             f'NovaSeq {pr(q[F]["reps|NovaSeq 30x"])} and G400RS 30x {pr(q[F]["reps|DNBSeq-G400RS 30x"])} in '
                             f'the replicates; DNBSeq-T7 {pr(q[F]["t7_original"])} original, '
                             f'{pr(q[F]["reps|DNBSeq-T7 30x"])} in the replicates; controls {pr(q[F]["controls"])}).'),
        ('Negative control', f'The baseline must change nothing: concordance {_f(G["G6"]["recall_vs_original|deep"]["concordance"])} '
                             f'with the original genotypes (G6) and {fmt_int(G["G9"]["calls"])} HLA calls identical (G9); '
                             f'{A["fidelity_controls"]} HiSeq X fidelity controls go through the same new processing.'),
        ('Blind spot', f'The G400RS 15x cases ({d["DNBseq-G400RS 15x"]["depth"]:.1f}x: depth cannot be added); a platform '
                       'effect that is not a depth effect, such as read length or chemistry (the DNBSeq-T7 cases lose only '
                       f'{_pct(1 - DS["t7_fraction"], 1)} of their reads). Covered by Q2 and Q3.'),
        ('Processing chain', 'T7_FRACTIONS → ARM_PLAN → EXTRACT_READS → SUBSAMPLE → CALL_GVCF → JOINT → GENOTYPE_QC → '
                             'CONCAT_QC → SWAP → RETEST_GWAS, RETEST_GENE; HLA_SUBSAMPLE → HLAHD → HLA_INPUTS → HLA_CHAIN → '
                             'HLA_MATRICES → RETEST_HLA → RETEST_OMNIBUS; CONCORDANCE, QC_REASONS, HLA_CONCORDANCE → '
                             'GATES_DS; DS_SUMMARY'),
        ('Output tables', CHECK_TABLES('depth') + '; `07.downsampling/<label>/{concordance,qc_reasons,hla}/`; '
                          '`07.downsampling/<label>/summary/ds.{exposure,qc_reasons}.tsv`'),
    ])
    cards['Q2'] = OrderedDict([
        ('Question', 'Do the four case platforms agree, and does any one carry the signal? If one platform\'s calls made it, '
                     'the signal would sit in that platform\'s cases and vanish without them (threat T2).'),
        ('Factor varied', 'Stratified: the cases by platform. Per cohort the pooled set, four strata (one platform\'s cases) and '
                          f'four leave-one-platform-out sets: {SU["platform"]} subsets, each with a separate null model.'),
        ('Held constant', 'Every control in every subset; the original genotypes; each signal\'s model (sex and 10 principal '
                          'components), GRM markers and threshold.'),
        ('Perturbed (N / I / F)', 'Cases by platform: ' + '; '.join(f'{PLAT_SHORT[p]} {t3(lambda c, p=p: fmt_int(cp[c][p]))}'
                                                                    for p in CASE_PLATFORMS) + f'; controls {ctrls3} in every subset.'),
        ('Comparator', 'The strata against each other (heterogeneity); the random-effects estimate against the pooled one; '
                       'each platform left out against z_exp = pooled β / SE without it; the difference β_k − β_(−k).'),
        ('Model and statistic', 'Per subset its own null model and step 2. Across strata: GLS with the shared-control '
                                'correlation ρ_kl = √(n_k n_l / ((n_k + n_0)(n_l + n_0))) (Lin & Sullivan 2009); Q (Cochran '
                                '1954) and I² (Higgins & Thompson 2002); REML with Hartung–Knapp (Hartung & Knapp 2001; '
                                'IntHout 2014) against the pooled estimate (Evangelou & Ioannidis 2013), at MAF ≥ 5 %; case-only multinomial '
                                'LRT (Begg & Zhang 1994). Genes: CMC burden β per subset.'),
        ('Readout', 'Q, its P read empirically against null SNPs of the same MAF bin and type, I²; the random-effects OR '
                    'against the pooled OR; z observed against z expected and P_diff.'),
        ('Decision rule', f'platform_Q: {CK.RULES["platform_Q"]} (the threshold is {N["q_threshold"]}); platform_RE: '
                          f'{CK.RULES["platform_RE"]}; lopo: {CK.RULES["lopo"]}. A flag in any one cohort.'),
        ('Positive control', 'None by design (no signal with a known platform effect); each stratum carries its power to '
                             'detect heterogeneity if that stratum had no effect. Post hoc, the strata flag the ADARB2 lead, '
                             'whose deletion the DNBSeq-T7 calls miss although the T7 reads carry it (deck CRAM check, 2026-11-06). The testability bar is not a '
                             'control: it sets aside readings where a large group (DNBSeq-T7) is left out; a small group '
                             'passes it and barely moves z.'),
        ('Negative control', f'{fmt_int(T["null_snps"])} random LD-pruned SNPs and indels by MAF bin: Q P < 0.05 in '
                             f'{t3(lambda c: _pct(N["platform_null"][c]))}; the analytic ρ against the correlation of z over the '
                             f'null SNPs ({_f(N["rho"]["analytic"][0])}–{_f(N["rho"]["analytic"][1])} against '
                             f'{_f(N["rho"]["null_snps"][0])}–{_f(N["rho"]["null_snps"][1])}, three cohorts); G5: the pooled '
                             're-fit reproduces the null model of the main analysis.'),
        ('Blind spot', 'DNBSeq-T7 holds most of the cases, so leaving them out is testable only for strong signals '
                       '(covered by Q3); genes (0–7 carriers per stratum: no per-stratum estimate); small strata often "not '
                       'testable"; an artifact confined to a small stratum barely moves z without it and is read by '
                       'P_diff (§5).'),
        ('Processing chain', 'MAKE_SUBSETS → FIT_NULL → STEP2_SUBSET (with PLATFORM_GENO) → PLATFORM_STATS → PLATFORM_REPORT; '
                             'PLATFORM_GENE, CARRIER_FIRTH, PLATFORM_OMNIBUS; GATE_G5'),
        ('Output tables', CHECK_TABLES('platform_Q', 'platform_RE', 'lopo') +
         '; `10.platform/<cohort>/platform.{strata,case_only,rho,null_calibration,null_Q}.tsv`'),
    ])
    cards['Q3'] = OrderedDict([
        ('Question', 'Does a second technology give the same answer in the same people? If a DNBSeq-T7 calling error made '
                     'the signal, the arrays, which type the DNBSeq-T7 cases and the HiSeq X controls, would disagree '
                     '(threat T3; also T4).'),
        ('Factor varied', 'Replaced: the WGS genotypes by SNP-array data of the same people, one marker per signal: the '
                          'lead typed, else imputed (R² ≥ 0.8), else a typed LD proxy (r² ≥ 0.8); DQB1:57:V imputed '
                          'with HIBAG\'s published Asian model, which never saw our WGS typing; gene drivers (full '
                          'mainland only): their WGS case carriers looked up on the imputed arrays. '
                          'check/array_crosscheck\'s comparison is re-run in each cohort with that cohort\'s '
                          'unrelated members, phenotype and covariates; G10: the full-mainland re-run reproduces it.'),
        ('Held constant', 'The people, the marker and the allele; the shared sites, the association model (sex and 10 '
                          'principal components) and the procedures of check/array_crosscheck (read only).'),
        ('Perturbed (N / I / F)', 'People with arrays ' + t3(lambda c: fmt_int(AR[c]['people'])) + ': DNBSeq-T7 cases '
                                  f'(ASA24v1) {t3(lambda c: fmt_int(AR[c]["cases"]))}, HiSeq X controls (JSA24v1) '
                                  f'{t3(lambda c: fmt_int(AR[c]["controls"]))}; the NovaSeq and G400RS cases have no array.'),
        ('Comparator', 'The WGS result in the same people, at the same marker and allele.'),
        ('Model and statistic', f'Typed on both sides: the array association against the WGS one. Imputed: the arrays were '
                                f'imputed with {IM["tool"]} (v{IM["version"]}, per the header of the imputed VCF) against '
                                f'{IM["panel"]}; that panel holds the controls\' own genomes, so only the cases\' imputed '
                                'genotypes are independent of WGS. Gene drivers: the WGS case carriers that the imputed arrays '
                                'confirm (dosage ≥ 0.5), counted once in full mainland.'),
        ('Readout', 'GWAS leads and DQB1:57:V, per cohort: z on WGS and on the array, consistent, not checkable, or '
                    'contradicted. Genes: case carriers confirmed of the WGS carriers, one cell for the three '
                    'cohorts (a description, no verdict).'),
        ('Decision rule', CK.RULES['array'] + '; contradicted = opposite direction or |z_WGS − z_array| > 3; the array P '
                                              'is not required. A flag in any one cohort. Genes: no rule (a '
                                              'missed rare carrier is expected of imputation).'),
        ('Positive control', 'Genome-wide sites with |z_WGS − z_array| > 4: ' + t3(lambda c: str(AR[c]['abs_dz_gt4'])) +
                             '; the five in full mainland are array-side calling errors in cases (check/array_crosscheck). '
                             'The sites with WGS P < 10⁻⁴ that show the same direction with array P < 0.01: ' +
                             t3(lambda c: f'{AR[c]["of_which_array_same_direction_p_lt_0.01"]} of '
                                          f'{AR[c]["sites_wgs_p_lt_1e-4"]}') +
                             f'; the z correlation is {t3(lambda c: _f(AR[c]["r_z_wgs_array"]))}.'),
        ('Negative control', 'Concordance ' + t3(lambda c: str(AR[c]['concordance_cases_T7_vs_ASA24v1'])) + ' in cases and ' +
                             t3(lambda c: str(AR[c]['concordance_controls_HiSeqX_vs_JSA24v1'])) + ' in controls over '
                             f'{fmt_int(AR[F]["shared_sites"])} shared sites. (λ ' +
                             t3(lambda c: f'{_f(AR[c]["lambda_wgs"])} / {_f(AR[c]["lambda_array"])}') + ', WGS / array, is '
                             'inflation both technologies share, not a negative control.)'),
        ('Blind spot', 'The NovaSeq and G400RS cases; imputed controls (they copy their own genomes); indels without a proxy '
                       '(the ADARB2 lead); rare carriers the imputation misses, and an acquired variant (JAK2 p.Val617Phe) '
                       'absent from the imputed arrays; an error shared by WGS and arrays: the arrays are confounded with '
                       'case status too, by another technology.'),
        ('Processing chain', 'ARRAY_COHORTS (cached; array_cohorts.py: check/array_crosscheck\'s compare, rare-direct '
                             'and imputed steps re-run per cohort from its inputs, gate G10; array_hibag.py: DQB1 with '
                             'HIBAG) → ARRAY_EVIDENCE (array_evidence.py: the verdict per signal and the gene carrier '
                             'counts)'),
        ('Output tables', CHECK_TABLES('array') + '; `11.array/array.{summary,gate}.json`'),
    ])
    cards['Q4'] = OrderedDict([
        ('Question', 'Does any one of the five control sources carry a signal? If one source made it, that source would stand '
                     'out, and leaving it out would break the signal (threat T5).'),
        ('Factor varied', 'Stratified: the controls by source (the info sheet\'s Project): every case against one source\'s '
                          f'controls and against all sources but one: {SU["sources"]} subsets per cohort, each with its own '
                          'null model.'),
        ('Held constant', 'Every case in every subset; the original genotypes; each signal\'s main model and threshold.'),
        ('Perturbed (N / I / F)', 'Controls by source: ' + '; '.join(
            f'{SRC_LABEL[s_]} {t3(lambda c, s_=s_: fmt_int(src[c][s_]))}' for s_ in SOURCES) + f'; cases {cases3} in every subset.'),
        ('Comparator', 'The sources against each other; each source left out against z_exp; one source against the rest.'),
        ('Model and statistic', 'As Q2 with cases and controls swapped: GLS with ρ_kl = √(m_k m_l / ((m_k + n)(m_l + n))), m a '
                                'source\'s controls, n the cases (Lin & Sullivan 2009); Q with its analytic χ² P, I²; genes '
                                'CMC and SKAT-O per subset (carriers per source descriptive); positions: the omnibus.'),
        ('Readout', 'Q P and I² across the five sources; z observed against z expected without each source; the genes\' P, OR '
                    'and carriers per subset.'),
        ('Decision rule', f'sources_Q: {CK.RULES["sources_Q"]}; sources_loo: {CK.RULES["sources_loo"]}. For a rare variant in a '
                          'small source subset the one-step β is unreliable: leaving a source out is read first.'),
        ('Positive control', 'None by design (no source with a known effect); every GWAS and HLA reading is testable '
                             'here, since all cases stay.'),
        ('Negative control', f'λ over the {fmt_int(T["null_snps"])} null SNPs in every subset, '
                             f'{_f(N["sources_lambda"][0])}–{_f(N["sources_lambda"][1])}; |z| > 3 in '
                             f'{100 * N["sources_share_z3"][0]:.2f}%–{100 * N["sources_share_z3"][1]:.2f}% of them '
                             f'({100 * N["share_z3_null"]:.2f}% expected).'),
        ('Blind spot', 'A difference shared by all five sources, such as their common HiSeq X platform (covered by Q3); genes: '
                       '0–3 carriers per source, no per-source estimate. Pooling controls from different projects is a known '
                       'source of false positives (Lee 2017; Guo 2018); that is what the check targets.'),
        ('Processing chain', 'MAKE_SUBSETS → MAKE_CTRL_SUBSETS → FIT_NULL → STEP2_SUBSET; PLATFORM_GENE; PLATFORM_OMNIBUS; '
                             'CONTROL_REPORT'),
        ('Output tables', CHECK_TABLES('sources_Q', 'sources_loo') + '; `12.controls/controls.{counts,strata,null_lambda,'
                                                                    'signals}.tsv`'),
    ])
    return cards


def CHECK_TABLES(*keys):
    """the registry's tables of these checks as code spans, entries with one file prefix merged"""
    groups = OrderedDict()
    for k in keys:
        for t in CK.CHECK[k].tables.split('; '):
            m = re.match(r'^(.*)\.\{(.+)\}(\.\w+)$', t) or re.match(r'^(.*)\.([^./{}]+)(\.\w+)$', t)
            pre, items, ext = (m.group(1), m.group(2).split(','), m.group(3)) if m else (t, [], '')
            g = groups.setdefault((pre, ext), [])
            g += [i for i in items if i not in g]
    return '; '.join(f'`{pre}.{{{",".join(it)}}}{ext}`' if len(it) > 1 else f'`{pre}.{it[0]}{ext}`' if it else f'`{pre}`'
                     for (pre, ext), it in groups.items())


def _methods(checks):
    return ', '.join(dict.fromkeys(x for c in checks for x in CK.CHECK[c].methods.split(', ')))


def design_markdown(N):
    F, A, M, R, C, G = 'full_mainland', N['arms'], N['matrix'], N['registry'], N['compute'], N['gates']
    DS, AR, FU, IM = N['downsampled'], N['array'], N['followup'], N['imputation']
    t3 = lambda f: ' / '.join(f(c) for c in COHORTS)  # noqa: E731
    ng = R['gates']
    L = ['# Design — `check/platform_depth`', '',
         'How the validation is designed: the threat model, the fidelity layer, four questions (Q1–Q4), each answered by one '
         'experiment in one template, their '
         'controls and coverage, the decision rules, and the units at every step. Written by `scripts/design_figures.py` '
         'from the registry (`scripts/checks.py`), the result tables (`results/`) and the run traces; do not edit by hand '
         '(`--check-doc` fails when it is out of date). Formulas: [METHODS.md](METHODS.md); processes and checks: '
         '[PIPELINE.md](PIPELINE.md); files: [OUTPUTS.md](OUTPUTS.md); results: [README](../README.md#results). Part II '
         'of the deck (`deck/`) draws the same figures and cards.', '',
         'Contents: [1 At a glance](#1-at-a-glance) · [2 Threat model](#2-threat-model) · '
         f'[3 Fidelity layer](#3-fidelity-layer-gates-g1g{ng}) · '
         '[4 Questions](#4-questions-one-experiment-each-one-template) · [5 Controls](#5-positive-and-negative-controls) · '
         '[6 Coverage](#6-coverage) · [7 Decision](#7-from-readouts-to-a-decision) · [8 Numbers flow](#8-numbers-flow) · '
         '[9 Process graph](#9-process-graph) · [10 Limitations](#10-limitations) · [11 Figures](#11-how-the-figures-are-made)', '']
    # 1
    L += ['## 1. At a glance', '', '![master workflow](design/master_workflow.png)', '',
          'Platform coincides with case status: every control was sequenced on HiSeq X (151-bp reads, '
          f'{N["depth"][CONTROL_PLATFORM]["depth"]:.1f}x), every case on one of four other platforms. A platform term would '
          'absorb case status, so each way an artifact could arise is tested by design instead. The main study and its '
          f'signals first pass the fidelity gates G1–G{ng}: the re-test machinery reproduces every main-analysis result before '
          'anything is compared, so a later change can only come from the factor an experiment varies. Four questions, '
          'Q1 depth, Q2 platform, Q3 technology and Q4 control source, are each answered by one experiment that varies '
          'one factor, in all three nested cohorts; each re-tests or reads every signal with its own model and threshold, '
          'and feeds '
          f'named checks: the columns of the evidence matrix ({M["signals"]} signals × {M["checks"]} checks × 3 cohorts, '
          f'{fmt_int(M["signals"] * M["checks"] * 3)} cells: {fmt_int(M["readings"])} defined readings, {fmt_int(M["applying"])} '
          'that apply). Fixed rules read the matrix; a flag never removes a signal; the group decides. "Consistent" means '
          'not flagged, not shown robust. The cohorts share most of their samples: agreement across them shows stability to '
          'sample selection, not replication.', '']
    # 2
    L += ['## 2. Threat model', '', '![threat model](design/threat_model.png)', '',
          'Each row is a way an artifact could reach the association; its prediction names what an experiment must see if '
          'the artifact made the signal. "Tests" = the experiment can detect it directly; "partly" = only part of the '
          'mechanism; which samples each experiment involves is shown in §6.', '']
    rows = []
    for t in threats(N):
        cov = COVERAGE[t['key']]
        rows.append([t['key'], t['difference'], t['evidence'], t['path'], t['prediction'],
                     ', '.join(EXP_NAME[e] if e == 'G' else e for e, v in cov.items() if v == 'tests'),
                     ', '.join(e for e, v in cov.items() if v == 'partly') or '—'])
    L += _md_table(['', 'Difference between sequencing groups', 'Evidence', 'Path to a false association',
                    'Testable prediction', 'Tests', 'Partly'], rows) + ['']
    # 3
    g1, g5 = [G['G1'][c] for c in COHORTS], [G['G5'][c] for c in COHORTS]
    g6, g6f = G['G6']['recall_vs_original|deep'], G['G6']['recall_vs_original|fidelity']
    res = {
        'G1': f'{t3(lambda c: fmt_int(G["G1"][c]["variants_ours"]))} variants; largest relative difference '
              f'{max(x["max_rel_diff"] for x in g1):g}',
        'G2': f'{len(G["G2"])} cohort × variant-set scans; mismatching rows {sum(len(x["mismatches"]) for x in G["G2"].values())}',
        'G3': f'{G["G3"]["identical"]} of {G["G3"]["files"]} files identical',
        'G4': f'{fmt_int(G["G4"]["single_rows"])} single-marker and {fmt_int(G["G4"]["omnibus_rows"])} omnibus rows; '
              f'largest relative difference {G["G4"]["worst"]:g}',
        'G5': f'variance ratio within {_sci(max(x["variance_ratio_rel_diff"] for x in g5))}, z within '
              f'{_sci(max(x["z_max_rel_diff"] for x in g5))}',
        'G6': f'concordance {_f(g6["concordance"])}: {g6["samples"]} down-sampled cases ({fmt_int(g6["both_called"])} '
              f'genotypes), {g6f["samples"]} HiSeq X controls; call rate {_f(g6["callrate_ref"], 5)} original, '
              f'{_f(g6["callrate_test"], 5)} re-processed',
        'G7': f'{A["swap_records"]} swap records; any other change aborts',
        'G8': f'{fmt_int(G["G8"]["genome_arms"])} genome-arms: median {_f(G["G8"]["ratio_median"])} '
              f'({_f(G["G8"]["ratio_min"])}–{_f(G["G8"]["ratio_max"])})',
        'G9': f'{fmt_int(G["G9"]["calls"])} calls; differing {G["G9"]["differing"]}',
        'G10': f'{G["G10"]["tables"]} tables; largest relative difference {_sci(G["G10"]["max_rel_diff"])}',
    }
    L += [f'## 3. Fidelity layer (gates G1–G{ng})', '',
          'Nothing is compared before the machinery reproduces the main-analysis results; the baseline (f = 1: every read, '
          'through the same new processing) reproduces the '
          'original genotypes and typing, and the array comparison re-run in full mainland reproduces check/array_crosscheck. '
          'A difference introduced by the re-test machinery would otherwise be read as an artifact (threat T0). Blind spot: '
          'a flaw shared by the original and the re-test (same software, same settings) passes the gates. '
          f'All {ng} gates pass{"" if G["pass"] else " — NOT: see results/_gates"}.', '']
    L += _md_table(['Gate', 'What is reproduced', 'Rule', 'Process', 'Record', 'Result'],
                   [[g.key, g.what, g.rule, f'`{g.process}`', f'`{g.record}`', res[g.key]] for g in CK.GATES]) + ['']
    # 4
    L += ['## 4. Questions, one experiment each, one template', '',
          'Every experiment is described by the same fields. "Perturbed" gives the narrow / intermediate / full mainland '
          'counts (N / I / F). Each feeds the evidence-matrix checks of the registry named in its heading; every check is '
          'read in all three cohorts.', '']
    cards = experiment_cards(N)
    for (k, name, checks, verb), (_, card) in zip(EXPERIMENTS, cards.items()):
        L += [f'### {k} · {name} — {verb}', '',
              f'Checks: {", ".join(f"`{c}`" for c in checks)} · METHODS {_methods(checks)}', '']
        L += _md_table(['Field', 'Content'], [[f'**{f}**', v] for f, v in card.items()]) + ['']
    # 5
    L += ['## 5. Positive and negative controls', '',
          'A positive control shows that the manipulation or the test could detect what it targets; a negative control shows '
          'that it changes nothing, or stays calibrated, where nothing should change.', '']
    L += _md_table(['Experiment', 'Positive control', 'Negative control'],
                   [[f'{k} · {n}', cards[k]['Positive control'], cards[k]['Negative control']] for k, n, *_ in EXPERIMENTS]) + ['']
    L += ['The leave-one-out rule compares z without a group with z expected = β_pooled / SE without the group. Given the '
          'pooled estimate, z without the group follows N(z_exp, sd = √(1 − (z_exp / z_pooled)²)): the winner\'s curse '
          'inflates z_pooled and z_exp alike and does not bias the comparison. The sd is '
          f'{_f(N["loo_sd"]["t7"][0], 2)}–{_f(N["loo_sd"]["t7"][1], 2)} without the DNBSeq-T7 cases and '
          f'{_f(N["loo_sd"]["small"][0], 2)}–{_f(N["loo_sd"]["small"][1], 2)} without a small group (GWAS and HLA markers); '
          'the testability bar excludes large groups (low z_exp) only, while a '
          'small group always passes it and barely moves z, so an excess confined to a small group is only described, by '
          'P_diff and, for a rare lead (whose one-step stratum β is unreliable), by its allele counts. On the null SNPs the '
          'model holds for every platform left out; for the control sources it is approximate (each subset has its own '
          'null model), so their chances are nominal. '
          '`14.evidence/evidence.loo_calibration.tsv` gives the flags observed against those expected by chance.', '']
    # 6
    L += ['## 6. Coverage', '', '![coverage](design/coverage.png)', '',
          'Each strip is every genome of full mainland on one scale, cases by platform and controls by source; color = what '
          'the experiment does to that group. No single experiment sees every sample: the DNBSeq-T7 and deep cases are '
          f'covered by at least two, the {DS["untouched"]} G400RS 15x cases by the platform strata alone.', '']
    head = ['Threat', 'Gates'] + [k for k, *_ in EXPERIMENTS]
    sym = {'tests': '●', 'partly': '◐'}
    L += _md_table(head, [[f'{t["key"]} {t["short"]}'] + [sym.get(COVERAGE[t['key']].get(c), '·') for c in ['G'] +
                                                          [k for k, *_ in EXPERIMENTS]] for t in threats(N)])
    L += ['', '● tests it directly · ◐ partly (only part of the mechanism) · · blind. How the blind spots are covered: a '
          'DNBSeq-T7 calling error (T3) is seen by the strata only where leaving the DNBSeq-T7 cases out is testable (they '
          'are most of the cases); the arrays (Q3) test it directly in the same people with a second technology. Depth is varied in every case platform but G400RS 15x ({DS["untouched"]} cases '
          'below the controls\' depth: depth cannot be added); a read-length or chemistry effect is not a depth effect and '
          'is left to the strata and the arrays. Imputed arrays cannot test controls (their panel is the controls\' own '
          'genomes); a control-side error shared by all five sources is not separable from population or panel '
          'differences. Indels without a proxy escape the arrays: the ADARB2 lead is seen by the strata only.', '']
    # 7
    L += ['## 7. From readouts to a decision', '',
          '1. **Reading.** Each signal × check × cohort is read by the check\'s fixed rule. ' +
          '; '.join(f'*{k}*: {v}' for k, v in CK.READINGS.items()) + '.',
          '2. **Cohort-wise rule.** Every check is read in each of the three cohorts and flags a signal when its rule holds '
          'in any one of them; the number of cohorts flagged is reported with the flag.',
          '3. **Flag.** The signal carries the flag in the matrix and in `14.evidence/evidence.signals.tsv`; a flag never '
          'removes a signal.',
          '4. **Look behind.** P_diff and the descriptive follow-up counts (`13.followup`) are read with every flag.',
          '5. **Decision.** Which flag keeps a locus out of the paper is the group\'s decision (no signal is excluded for '
          'now). The rules were fixed before the final readings and are to be confirmed by the group. Three changes came '
          f'after review: the random-effects rule was restricted, like Q, to MAF ≥ {_pct(float(R["q_min_maf"]), 0)}; '
          f'{A["replicates"]} replicates for every platform was the group\'s choice after the first run; and every check '
          'flags a signal when its rule holds in any one cohort, not two of three (2026-11-04: the cohorts are nested, so '
          'a second cohort adds little, and a flag only prompts a closer look).', '',
          f'Thresholds (registry): leave-one-out |z expected| ≥ {R["loo_z_expected"]}, |z observed| < {R["loo_z_observed"]}; '
          f'depth Δ|z| > {R["depth_dz"]}; heterogeneity {R["q_alpha"]} / the number of SAIGE signals at MAF ≥ '
          f'{R["q_min_maf"]} (= {N["q_threshold"]}).', '']
    L += _md_table(['Check', 'Experiment', 'Rule (`evidence.rules.tsv`)'],
                   [[f'`{c.key}` {c.title}', next((k for k, _, ch, _ in EXPERIMENTS if c.key in ch), 'main analysis'), c.rule]
                    for c in CK.CHECKS]) + ['']
    # 8
    S, SU, T = N['signals'], N['subsets'], N['targets']
    hla = {c: S['hla_alleles'][c] + S['hla_residues'][c] for c in COHORTS}
    L += ['## 8. Numbers flow', '', '![numbers flow](design/numbers_flow.png)', '']
    L += _md_table(['Step', 'Units'], [
        ['Cohorts (N / I / F)', f'genomes {t3(lambda c: fmt_int(N["genomes"][c]))}; cases {t3(lambda c: fmt_int(N["cases"][c]))}; '
                                f'controls {t3(lambda c: fmt_int(N["controls"][c]))}'],
        ['Q1 genomes', f'{DS["n"]} down-sampled cases: {DS["deep"]} deep cases (the same in every cohort) and {DS["t7"]} '
                       f'DNBSeq-T7 cases (f = {DS["t7_fraction"]}); per cohort {t3(lambda c: fmt_int(DS["by_cohort"][c]["all"]))}; '
                       f'{A["fidelity_controls"]} HiSeq X fidelity controls'],
        ['Q1 arms', f'{A["arms"]} arms (baseline and {A["replicates"]} replicates): {fmt_int(A["deep_genome_arms"])} case '
                    f'genome-arms; {fmt_int(A["called_genome_arms"])} re-called with the fidelity controls'],
        ['Q1 joint call', f'{fmt_int(A["joint_gvcfs"])} GVCFs ({A["original_gvcfs"]} original) at {fmt_int(T["fixed_sites"])} '
                          f'fixed sites ({T["by_kind"]["gwas_signal"]} signal, {fmt_int(T["by_kind"]["gwas_context"])} context, '
                          f'{T["by_kind"]["gene"]} gene) in {T["intervals"]} calling intervals ({T["interval_kb"]:.0f} kb)'],
        ['Q1 HLA and swaps', f'{fmt_int(A["hla_typings"])} HLA re-typings; {A["swap_records"]} swap records'],
        ['Q1 re-tests', f'{A["retest_arms"]} arms × 3 cohorts: {rng(S["gwas_retest_variants"].values())} variants, '
                        f'{S["gene_sets"]} gene sets, {rng(hla.values())} HLA markers, {rng(S["positions"].values())} positions'],
        ['Q2, Q4 subsets', f'{SU["per_cohort"]} per cohort ({SU["platform"]} platform, {SU["sources"]} control-source): '
                           f'{SU["null_models"]} null models'],
        ['Q2, Q4 tests', f'{SU["step2_runs"]} step-2 runs of {rng([x for v in SU["markers_per_run"].values() for x in v])} '
                         f'markers ({fmt_int(T["null_snps"])} of them null SNPs and indels): {fmt_int(SU["step2_tests"])} tests; '
                         f'{SU["gene_runs"]} gene runs; {SU["omnibus_runs"]} omnibus runs'],
        ['Q3', f'people with arrays {t3(lambda c: fmt_int(AR[c]["people"]))} (cases {t3(lambda c: fmt_int(AR[c]["cases"]))}, '
               f'controls {t3(lambda c: fmt_int(AR[c]["controls"]))}); {fmt_int(AR[F]["shared_sites"])} shared sites; '
               f'{AR[F]["signals"]} signals per cohort; imputed with {IM["tool"]} v{IM["version"]}'],
        ['Follow-up counts', f'{t3(lambda c: str(FU[c]["gene_variants"]))} gene variants × {t3(lambda c: fmt_int(FU[c]["genomes"]))} '
                             f'genomes re-filtered; {S["loci"]} leads by group (descriptive, not a check)'],
        ['Signals', f'{S["loci"]} GWAS loci ({S["gwas_variants"]} signal variants: {S["leads"]} leads, {S["alt_leads"]} '
                    f'alternative leads, {S["cs_members"]} credible-set members); {len(S["genes"])} Tier 1–2 genes in '
                    f'{S["gene_sets"]} gene sets ({T["gene_variants"]} tested variants); HLA '
                    f'{rng(S["hla_alleles"].values())} alleles and {rng(S["hla_residues"].values())} residues, '
                    f'{rng(S["positions"].values())} positions'],
        ['Matrix', f'{M["signals"]} signals × {M["checks"]} checks × 3 cohorts = {fmt_int(M["signals"] * M["checks"] * 3)} '
                   f'cells: {fmt_int(M["readings"])} defined readings, {fmt_int(M["applying"])} that apply (and '
                   f'{M["discovery_readings"]} discovery readings); {M["flag_cells"]} flag readings in '
                   f'{len(M["flag_signals"])} signals ({", ".join(M["flag_signals"])}); none removed'],
        ['Compute', f'behind these results {fmt_int(C["tasks"])} tasks of {C["processes"]} processes, '
                    f'{C["core_h"]:,.1f} core-hours allocated ({C["cpu_h"]:,.1f} CPU-hours used); by stage ' +
                    ', '.join(f'{k} {v:,.1f}' for k, v in C["core_h_by_stage"].items()) +
                    f'; all {C["runs"]} runs with the pilots: {fmt_int(C["executed_tasks"])} tasks, '
                    f'{C["executed_core_h"]:,.1f} core-hours.'],
    ]) + ['', 'Core-hours: each task\'s realtime × the cores `nextflow.config` gives its process (`--rsc c=`); "used" = '
          'realtime × %cpu. Tasks behind these results: those the latest full run and the latest report run used, '
          'each counted once at the run that executed it.', '',
          f'Replicates: {A["replicates"]} for every platform, our choice for uniformity. The pilot measured the cost '
          'of every down-sampling process per task; its cost rule (the largest R in 3–10 within 1,000 core-hours and 48 h '
          f'at 200 cores) gave R = {N["pilot"]["R_rule"]} for the deep cases alone.', '']
    # 9
    L += ['## 9. Process graph', '',
          f'The {R["processes"]} documented processes by stage (' + ', '.join(f'{k} {v}' for k, v in R['stages'].items()) +
          '). Edges: the channels of `platform_depth.nf` (parsed by `channel_edges()`); the report stage reads the published '
          'folders named by the registry. Green: a process that writes a gate record. The sample table and the targets feed '
          'most processes and are drawn to each stage. Every edge in one graph: [design/process_dag.svg](design/process_dag.svg).', '',
          '![process graph, part 1](design/process_dag_part1.png)', '', '![process graph, part 2](design/process_dag_part2.png)', '']
    # 10
    L += ['## 10. Limitations', '',
          '- Platform and case status coincide: no experiment can remove the confounding; each tests one consequence of it.',
          '- A flaw shared by the original analysis and the re-test passes the gates.',
          '- Depth is varied in every case platform but G400RS 15x (below the controls\' depth); the DNBSeq-T7 cases lose '
          f'only {_pct(1 - DS["t7_fraction"], 1)} of their reads, so a read-length or chemistry effect is left to the strata '
          'and the arrays.',
          '- A DNBSeq-T7 calling error is tested directly only by the arrays, which exist for the DNBSeq-T7 cases and the '
          'HiSeq X controls only, and are themselves confounded with case status (by another technology).',
          '- Small strata and sources are often "not testable"; genes are read by leaving groups out only (too few carriers '
          'per group for an estimate).',
          '- Leaving out a small group barely moves z: an artifact confined to a small group is read by P_diff, not by the '
          'leave-one-out rule.',
          '- The follow-up counts (carriers by group, possible carriers set missing by the genotype filters) describe and do '
          'not test. The matrix rules were fixed before the final readings and are to be confirmed by the group.', '']
    # 11
    L += ['## 11. How the figures are made', '',
          '```bash', 'cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/platform_depth',
          'PYTHONDONTWRITEBYTECODE=1 /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -B scripts/design_figures.py',
          'PYTHONDONTWRITEBYTECODE=1 /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -B scripts/design_figures.py --check-doc',
          '```', '',
          'It writes `docs/design/*.png` (300 dpi), the graphviz sources `docs/design/*.dot` (render: `FONTCONFIG_FILE='
          'deck/fonts.conf dot -Gdpi=300 -Tpng -o x.png x.dot`), `docs/design/figures.json` (natural size and smallest font '
          'of each figure: the deck places them at that size, so every text stays at least 9 pt) and this file. Read only: '
          '`scripts/checks.py`, `platform_depth.nf`, `nextflow.config`, `results/`, `deck/references.json`. '
          '`deck/collect_stats.py` imports `numbers()`, so the deck prints the same values.', '']
    figs = json.loads((OUT / 'figures.json').read_text())
    L += _md_table(['Figure', 'Tool', 'Size (in)', 'Smallest text (pt)'],
                   [[f'`{k}.png`', v['tool'], f'{v["width_in"]:.1f} × {v["height_in"]:.1f}', v['min_font_pt']]
                    for k, v in figs.items()]) + ['']
    return '\n'.join(L)


if __name__ == '__main__':
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('--numbers', action='store_true', help='print numbers() as JSON and exit')
    ap_.add_argument('--edges', action='store_true', help='print the process edges and exit')
    ap_.add_argument('--check-doc', action='store_true', help='exit 1 if docs/DESIGN.md is not what this script writes')
    a = ap_.parse_args()
    if a.numbers:
        print(json.dumps(numbers(), indent=1, default=str))
    elif a.edges:
        E = channel_edges()
        for e in E:
            print(*e, sep='\t')
        print('folder reads:', folder_reads())
        print('registry mentions not among the channel edges:', sorted(registry_mentions() - set(E)))
    elif a.check_doc:
        ok = (ROOT / 'docs' / 'DESIGN.md').read_text() == design_markdown(numbers()) + '\n'
        print('docs/DESIGN.md', 'is up to date' if ok else 'differs from what scripts/design_figures.py writes')
        sys.exit(0 if ok else 1)
    else:
        N = numbers()
        threat_model(N)
        coverage(N)
        master_workflow(N)
        numbers_flow(N)
        for part in ('full', 'A', 'B'):
            process_dag(part)
        (OUT / 'figures.json').write_text(json.dumps(MANIFEST, indent=1) + '\n')
        (ROOT / 'docs' / 'DESIGN.md').write_text(design_markdown(N) + '\n')
        missing = sorted(registry_mentions() - set(channel_edges()))
        if missing:
            sys.exit(f'the registry names upstream processes the workflow does not connect: {missing}')
        print('wrote', ', '.join(f'docs/design/{k}.png' for k in MANIFEST), 'and docs/DESIGN.md')
