#!/bin/bash
# Pre-flight and post-run checks for gene_signals. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree and the upstream inputs and must pass before
# a run. Sections 6-9 read results/ (or $OUT) and are skipped until it exists;
# after a run they must pass too. Section 8 compares against anchors taken from
# the assoc_gene run of 2026-10 and runs only while its gene_scan_all.tsv is
# unchanged (md5 below).
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/gene_signals
cd "$C"
OUT=${OUT:-results}
ANCHOR_MD5=36f2aa8462b6835799236412ffc505e5
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
source activate dsl2 2>/dev/null
out=$(nextflow lint gene_signals.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
# the pipeline's own interpreter, first on PATH for every python3 below
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" gene_signals.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import xlsxwriter, openpyxl, pandas" 2>/dev/null
ck "params.python ($PYBIN) imports xlsxwriter, openpyxl, pandas" $?

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' gene_signals.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' gene_signals.nf; ck "no task.cpus (ext.threads only)" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E '^\s*clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M"; bad=1; }
done
ck "every clusterOptions m= is a multiple of 4571" $bad
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*(python3|Rscript|source activate)\s' gene_signals.nf; ck "no bare python3 / Rscript / source activate in a process" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*import\s' gene_signals.nf; ck "no import in the .nf (strict syntax)" $([ $? -eq 1 ] && echo 0 || echo 1)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('gene_signals.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
cfg = pathlib.Path('nextflow.config').read_text()
bad = 0
for p, b in bodies.items():
    if p in ('WRITE_RUN_MANIFEST', 'OPENTARGETS'):        # local: provenance / needs the internet
        ok = "executor 'local'" in b
    else:
        ok = "executor 'slurm'" in b and "queue 'gr10478b'" in b and re.search(r"\n\s*time '", b)
        if not re.search(rf"withName: '[^']*\b{p}\b", cfg):
            print(f'    ✗ {p}: no resources in nextflow.config'); bad += 1
    if not ok:
        print(f'    ✗ {p}: executor / queue / time'); bad += 1
    for m in re.finditer(r'publishDir\s+"([^"]*)"', b):
        for var in re.findall(r'\$\{([^}]+)\}', m.group(1)):
            if not var.strip().startswith('params.'):
                print(f'    ✗ {p}: publishDir without a closure references {var}'); bad += 1
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "executor/queue/time and resources per process; publishDir closures" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('gene_signals.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
PAIRS = {'SELECT_GENES': ('select_genes.py', None), 'GENOTYPE_STATS': ('genotype_stats.py', None),
         'RV_QC_STATS': ('rv_qc_stats.py', None), 'VARIANT_SETS': ('variant_sets.py', None),
         'MAP_RECOUNT': ('map_recount.py', None),
         'ANNOTATE_FUNCTION': ('annotate_function.py', None), 'REF_ANNOVAR': ('annovar_annotations.py', None),
         'REF_TOMMO': ('ref_frequencies.py', 'tommo'), 'REF_1KG': ('ref_frequencies.py', 'kg'),
         'OPENTARGETS': ('opentargets.py', None), 'ASSEMBLE': ('assemble_table.py', None),
         'WRITE_TABLE': ('write_table.py', None)}
bad = 0
# LOVO calls lovo.py twice, once per subcommand: each command's flags against that subcommand's --help
for sub in ('sets', 'collect'):
    m = re.search(r'\$\{params\.python\} \$\{script\} ' + sub + r' ((?:[^\n]*\\\n)*[^\n]*)', bodies.get('LOVO', ''))
    if not m:
        print(f'  ✗ LOVO: no lovo.py {sub} call'); bad += 1; continue
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    h = subprocess.run([sys.executable, 'scripts/lovo.py', sub, '--help'], capture_output=True, text=True)
    known = set(re.findall(r'--([a-z][a-z0-9-]*)', h.stdout))
    if h.returncode != 0 or passed - known:
        print(f'  ✗ LOVO: lovo.py {sub} does not accept --' + ', --'.join(sorted(passed - known))); bad += 1
    if known - passed - {'help'}:
        print(f'  ✗ LOVO: lovo.py {sub} flags not passed: --' + ', --'.join(sorted(known - passed - {'help'}))); bad += 1
for proc, body in bodies.items():
    if '${params.python} ${script}' in body and proc not in PAIRS and proc != 'LOVO':
        print(f'  ✗ {proc} runs a script that verify.sh does not know'); bad += 1
for proc, (sc, sub) in sorted(PAIRS.items()):
    m = re.search(r'\$\{params\.python\} \$\{script\}(.*?)"""', bodies.get(proc, ''), re.S)
    if not m:
        print(f'  ✗ {proc}: no ${{params.python}} ${{script}} line'); bad += 1; continue
    # flags written in the command, and those built in the workflow (ASSEMBLE's QC arguments)
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    if proc == 'VARIANT_SETS':                         # the --map-dir argument is built in the workflow
        passed |= set(re.findall(r"'--([a-z][a-z0-9-]*) '", re.search(r'def mapArgs = (.*?)\n', nf).group(1)))
    cmd = [sys.executable, f'scripts/{sc}'] + ([sub] if sub else []) + ['--help']
    h = subprocess.run(cmd, capture_output=True, text=True)
    if h.returncode != 0:
        print(f'  ✗ {proc}: {sc} --help fails'); bad += 1; continue
    known = set(re.findall(r'--([a-z][a-z0-9-]*)', h.stdout))
    miss = passed - known
    if miss:
        print(f'  ✗ {proc}: {sc} does not accept --' + ', --'.join(sorted(miss))); bad += 1
    req = set(re.findall(r'add_argument\(\'--([a-z][a-z0-9-]*)\', [^)]*required=True',
                         pathlib.Path(f'scripts/{sc}').read_text()))
    if sub is None and req - passed:
        print(f'  ✗ {proc}: {sc} requires --' + ', --'.join(sorted(req - passed))); bad += 1
print(f'    {len(PAIRS)} process/script pairs checked')
sys.exit(1 if bad else 0)
PY
ck "flags accepted, required flags passed" $?

echo "═══ 4. upstream inputs ═══"
P=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6
U=$P/analysis/assoc_gene/results
bad=0
for f in _run_info/run_manifest.json _comparison/robust_genes/tiers.tsv; do [ -s "$U/$f" ] || { echo "    ✗ $U/$f"; bad=1; }; done
for t in gene_scan_all gene_index_all variant_stats_all set_stats_all scan_qc_all denominators cohort_samples multi_carriers_all; do
  [ -s "$U/_comparison/tables/$t.tsv" ] || { echo "    ✗ $t.tsv"; bad=1; }; done
for c in narrow_mainland intermediate_mainland full_mainland; do
  for f in anno_counts.tsv map_qc.tsv map.moderate_high.tsv.gz map.low_moderate_high.tsv.gz gene_index.moderate_high.tsv \
           gene_index.low_moderate_high.tsv denominator.moderate_high.tsv denominator.low_moderate_high.tsv; do
    [ -s "$U/$c/01.map/$f" ] || { echo "    ✗ $U/$c/01.map/$f"; bad=1; }; done
done
ck "assoc_gene manifest, tiers, comparison tables and per-cohort gene maps present" $bad
python3 - "$U" <<'PY'
import json, sys
m = json.load(open(f'{sys.argv[1]}/_run_info/run_manifest.json'))
ok = m['component'] == 'assoc_gene' and set(m['methods']) == {'cmc', 'skato'} and \
     set(m['strata']) == {'moderate_high', 'low_moderate_high'}
print(f"    upstream: cohorts {m['cohorts']}, minAC {m['minac']}, covariates {m['covariates'].count(',') + 1} terms")
sys.exit(0 if ok else 1)
PY
ck "upstream is assoc_gene with cmc + skato and the two variant sets" $?
bad=0
for c in narrow_mainland intermediate_mainland full_mainland; do
  T=$P/tuning.rv/results/$c
  for f in 00.prep/pheno_rvt.tsv 00.prep/covar_rvt.tsv 00.prep/pheno_coding.json 01.map/variants.moderate_high.txt \
           01.map/variants.low_moderate_high.txt; do
    [ -s "$U/$c/$f" ] || { echo "    ✗ $U/$c/$f"; bad=1; }; done
  for f in 02.callset_filter/filtered.bed 02.callset_filter/filtered.fam 00.precheck/02_variant_qc/variants_depthdiff.exclude.txt \
           00.precheck/02_variant_qc/depthdiff_stats.tsv 00.precheck/02_variant_qc/variant_metrics.txt.gz \
           00.precheck/metrics/sample_metrics.txt.gz 00.precheck/01_sample_qc/detect_outlier_samples.log \
           00.precheck/01_sample_qc/sample_outliers.remove.txt 05.qc_collect/minac_recommendation.tsv \
           03.qc_metrics/minac1/sample_metrics.txt.gz 03.qc_metrics/minac2/sample_metrics.txt.gz; do
    [ -s "$T/$f" ] || { echo "    ✗ $T/$f"; bad=1; }; done
  for s in .maf_ge_threshold.variants.txt .maf_ref.afreq; do
    n=$(ls -1 $P/wgs.auto.par/results/12_model_inputs/$c/genotype/fixed_model/*$s 2>/dev/null | wc -l)
    [ "$n" -eq 1 ] || { echo "    ✗ $c: $n *$s"; bad=1; }; done
  [ -s "$P/wgs.auto.par/results/12_model_inputs/$c/covariates/bbj_mainland_pc.sex.tsv" ] || { echo "    ✗ $c covariates"; bad=1; }
done
[ -s "$P/tuning.rv/results/_comparison/all_cohorts_qc_stats.tsv" ] || { echo "    ✗ tuning sweep"; bad=1; }
[ -s "$P/tmp/cram.v6/cram.v6.summary.csv" ] || { echo "    ✗ sample sheet"; bad=1; }
for f in variant_annot.py gene_annotate.R gene_utils.R; do [ -s "$P/analysis/_shared/scripts/$f" ] || { echo "    ✗ shared $f"; bad=1; }; done
[ -s "$P/analysis/assoc_gene/scripts/build_gene_map.py" ] || { echo "    ✗ assoc_gene build_gene_map.py"; bad=1; }
for f in $(grep -oE "params\.(plink2_test|rvtest_bin)\s*=\s*'[^']+'" gene_signals.nf | grep -oE "'[^']+'$" | tr -d "'"); do
  [ -e "$f" ] || { echo "    ✗ $f"; bad=1; }; done
ck "tuning.rv QC files, rare/common split, covariates, sample sheet, shared scripts, test inputs, rvtest present" $bad

echo "═══ 5. docs ═══"
bad=0; for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md; do
  [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
ck "README, METHODS, OUTPUTS, STUDY_NOTES present" $bad

if [ ! -d "$OUT" ] || [ ! -f "$OUT/checks.json" ]; then
  echo "═══ post-run sections skipped: no $OUT/checks.json ═══"
  if [ $fail -eq 0 ]; then echo "═══ PRE-FLIGHT PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
  exit $fail
fi

echo "═══ 6. DAG cardinality (trace.txt) ═══"
OUT="$OUT" python3 - <<'PY'
import csv, json, os, sys, collections
out = os.environ['OUT']
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
n_coh = len(m['cohorts'])
n_set = len(m['strata'])
EXPECT = {'SELECT_GENES': 1, 'GENOTYPE_STATS': n_coh, 'RV_QC_STATS': n_coh, 'MAP_RECOUNT': n_coh, 'VARIANT_SETS': 1,
          'LOVO': n_coh * n_set,
          'ANNOTATE_FUNCTION': 1, 'REF_ANNOVAR': 1, 'REF_TOMMO': 1, 'REF_1KG': 1, 'OPENTARGETS': 1, 'ASSEMBLE': 1,
          'WRITE_TABLE': 1, 'WRITE_RUN_MANIFEST': 1}
got, bad = collections.Counter(), 0
for r in csv.DictReader(open(f'{out}/_run_info/trace.txt'), delimiter='\t'):
    p = r['name'].split(' ')[0]
    if r['status'] in ('COMPLETED', 'CACHED'):
        got[p] += 1
    else:
        print(f"    ✗ {r['name']}: {r['status']}"); bad += 1
for p, n in EXPECT.items():
    if got[p] != n:
        print(f'    ✗ {p}: {got[p]} task(s), expected {n}'); bad += 1
print(f'    {sum(got.values())} tasks')
sys.exit(1 if bad else 0)
PY
ck "every task COMPLETED/CACHED, cardinality as expected" $?

echo "═══ 7. checks.json and run manifest ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, sys
out = os.environ['OUT']
c = json.load(open(f'{out}/checks.json'))
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
s, g, v = c['selection'], c['genes'], c['variants']
n_coh, n_set = len(m['cohorts']), len(m['strata'])
ok = {'all_passed': c['all_passed'],
      'gene rows = genes x sets x cohorts': g['rows'] == s['n_genes'] * n_set * n_coh,
      'variant rows = gene-variant pairs x cohorts': v['rows'] == s['n_gene_variant_rows'] * n_coh,
      'minAC as upstream': all(c['qc'][k]['minac'] == m['upstream_minac'][k] for k in m['cohorts']),
      'no variant kept only by protection': all(c['qc'][k]['gap_protected_rescued'] == 0 for k in m['cohorts']),
      'absent variants explained': set(v['callset_absent']) <= {'no_carrier', 'common', 'gap_filter'},
      'QC flow per cohort': set(c['qc']) == set(m['cohorts']),
      'LOVO: a baseline equal to the published result per tested gene x set x cohort x statistic':
          c['lovo']['baselines_equal_to_published'] == 2 * g['tested'],
      'LOVO: every gene has a primary variant set': set(c['lovo']['primary_sets']) == set(s['genes_in_order'])}
for k, x in ok.items():
    if not x:
        print(f'    ✗ {k}')
print(f"    {s['n_genes']} genes, {s['n_variants']} variants; coverage {c['coverage']}")
sys.exit(0 if all(ok.values()) else 1)
PY
ck "all_passed; rows partition; minAC; reasons" $?

echo "═══ 8. anchors (assoc_gene run of 2026-10; tuning.rv run of 2026-10-08) ═══"
md5=$(md5sum "$U/_comparison/tables/gene_scan_all.tsv" | cut -d' ' -f1)
if [ "$md5" != "$ANCHOR_MD5" ]; then
  echo "    (upstream gene_scan_all.tsv changed: anchors do not apply)"
else
OUT="$OUT" python3 - <<'PY'
import glob, json, os, sys
import numpy as np, pandas as pd
out = os.environ['OUT']
c = json.load(open(f'{out}/checks.json'))
s = c['selection']
want = {
    'genes in order': (s['genes_in_order'], ['STBD1', 'OSBP', 'RPH3A', 'JAK2', 'PIK3C2G', 'CD248', 'SCAMP3', 'NKD2',
                                              'PES1', 'SCFD1', 'P2RY2', 'FAM217B', 'TBC1D14', 'FAS', 'EFCAB6']),
    'classes': (s['class_counts'], {'Tier 1': 2, 'Tier 2': 2, 'Tier 3': 7, 'B': 3, 'S': 1}),
    'variants': ((s['n_variants'], s['impact_counts']), (141, {'MODERATE': 75, 'LOW': 58, 'HIGH': 8})),
    'absent': (c['variants']['callset_absent'], {'no_carrier': 12, 'common': 7}),
    # sensitivity family (>= 3 variants): OSBP and P2RY2 leave the moderate/high set, NKD2 rises, EFCAB6 enters
    'sensitivity tiers': (c['sensitivity']['tiers'],
                          {'CD248|low_moderate_high': 3, 'EFCAB6|low_moderate_high': 3, 'JAK2|moderate_high': 2,
                           'NKD2|low_moderate_high': 3, 'NKD2|moderate_high': 2, 'PES1|low_moderate_high': 3,
                           'PES1|moderate_high': 3, 'PIK3C2G|low_moderate_high': 3, 'PIK3C2G|moderate_high': 3,
                           'RPH3A|low_moderate_high': 2, 'RPH3A|moderate_high': 2, 'SCAMP3|low_moderate_high': 3,
                           'SCFD1|moderate_high': 3, 'STBD1|moderate_high': 1}),
}
bad = 0
for k, (got, exp) in want.items():
    if got != exp:
        print(f'    ✗ {k}: {got} != {exp}'); bad += 1
g = pd.read_csv(f'{out}/gene_signals.genes.long.tsv', sep='\t')
x = g[(g.gene_symbol == 'STBD1') & (g.stratum == 'moderate_high')].set_index('cohort')
exp = {'narrow_mainland': (1.31634e-05, 2.13224e-05, 9, 4), 'intermediate_mainland': (7.1974e-06, 1.20885e-05, 10, 7),
       'full_mainland': (7.22605e-08, 1.41108e-07, 10, 7)}
for coh, (pc, ps, ca, ct) in exp.items():
    r = x.loc[coh]
    if not (np.isclose(r.cmc_p, pc, rtol=1e-5) and np.isclose(r.skato_p, ps, rtol=1e-5)
            and (r.carriers_case, r.carriers_ctrl) == (ca, ct)):
        print(f'    ✗ STBD1 moderate_high {coh}'); bad += 1
# the rare-variant QC table (as published in the QC summary, 2026-10-08)
qc = {json.load(open(f))['cohort']: json.load(open(f)) for f in glob.glob(f'{out}/05.rv_qc/rv_qc.*.json')}
exp_qc = {'narrow_mainland': (2168, 419, 2165, 16103343, 16058030, 15939305, 7431380),
          'intermediate_mainland': (2480, 429, 2477, 17484447, 17440277, 17394531, 8050528),
          'full_mainland': (3071, 439, 3068, 19665336, 19622878, 19610994, 9161666)}
for coh, e in exp_qc.items():
    q, vq = qc[coh], qc[coh]['variant_qc']
    got = (q['n_in']['all'], q['n_in']['case'], q['n_kept']['all'], vq['n_in'], vq['n_in'] - vq['n_mono_after_sample_qc'],
           q['sweep']['1'], q['sweep']['2'])
    if got != e:
        print(f'    ✗ QC {coh}: {got} != {e}'); bad += 1
vs = json.load(open(f'{out}/05.rv_qc/variant_sets.json'))
if (vs['narrow_mainland']['analysed'], vs['narrow_mainland']['assigned_to_gene']) != (7431380, 145225):
    print('    ✗ variant sets (narrow)'); bad += 1
# leave-one-variant-out: the primary variant set of each gene, and the variant driving STBD1 in full (M/H)
lo = c['lovo']
mh, lmh = ['moderate_high'], ['low_moderate_high']
exp_prim = {**{g: mh for g in ('STBD1', 'OSBP', 'JAK2', 'NKD2', 'SCFD1', 'P2RY2', 'FAM217B')},
            **{g: lmh for g in ('CD248', 'SCAMP3', 'TBC1D14', 'FAS', 'EFCAB6')},
            **{g: sorted(mh + lmh) for g in ('RPH3A', 'PIK3C2G', 'PES1')}}
if lo['primary_sets'] != exp_prim:
    print(f"    ✗ LOVO primary sets {lo['primary_sets']}"); bad += 1
for meth in ('cmc', 'skato'):
    if lo['top'].get(f'STBD1|moderate_high|full_mainland|{meth}') != 'chr4:76309559:G:T':
        print(f'    ✗ LOVO STBD1 full M/H {meth} top contributor'); bad += 1
L = pd.read_csv(f'{out}/gene_signals.lovo.tsv', sep='\t')
x = L[(L.set_name == 'STBD1') & (L.stratum == 'moderate_high') & (L.cohort == 'full_mainland') & (L.method == 'cmc')]
if not np.isclose(x[x.dropped == 'chr4:76309559:G:T'].p.iloc[0], 0.000577461, rtol=1e-5):
    print('    ✗ LOVO STBD1 full M/H CMC P without chr4:76309559:G:T'); bad += 1
print('    gene order, classes, absent-variant reasons, sensitivity tiers, STBD1 burden and P, QC table, variant '
      'sets, LOVO primary sets and STBD1 driver compared')
sys.exit(1 if bad else 0)
PY
ck "anchors reproduced" $?
fi

echo "═══ 9. workbook ═══"
OUT="$OUT" python3 - <<'PY'
import glob, json, os, re, sys, warnings, zipfile
import pandas as pd
from openpyxl import load_workbook
out = os.environ['OUT']
with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    wb = load_workbook(f'{out}/gene_signals.xlsx', rich_text=True)
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
want = ['Rare-variant QC', 'Variant sets', 'Consequence terms', 'Gene-based scans', 'Gene-based signals',
        'Gene-based variants', 'Leave-one-variant-out', 'Gene-based Open Targets',
        'Multi-site carriers']   # the order of the workbook
if not m.get('opentargets'):
    want.remove('Gene-based Open Targets')
bad = 0
if wb.sheetnames != want:
    print(f'    ✗ sheets {wb.sheetnames}'); bad += 1
z = zipfile.ZipFile(f'{out}/gene_signals.xlsx')
for name in z.namelist():
    if name.startswith('xl/worksheets/sheet') and 'inlineStr' in z.read(name).decode():
        print(f'    ✗ inline strings in {name}'); bad += 1
    if name.startswith('xl/worksheets/sheet'):
        for pane in re.findall(r'<pane [^>]*>', z.read(name).decode()):
            if 'xSplit=' not in pane and 'ySplit=' not in pane:
                print(f'    ✗ frozen pane without a split in {name} (Excel repair prompt)'); bad += 1
ids = set()
for f in glob.glob(f'{out}/01.genotypes/samples.*.tsv'):
    ids |= set(pd.read_csv(f, sep='\t', dtype=str).IID)
# sample IDs and project names belong to sheet Multi-site carriers only
text = ' '.join(str(c.value) for sh in wb.worksheets if sh.title != 'Multi-site carriers' for row in sh.iter_rows()
                for c in row if c.value is not None)
ms_text = ' '.join(str(c.value) for row in wb['Multi-site carriers'].iter_rows() for c in row if c.value is not None)
for tok in ('_mainland', 'moderate_high', 'set_name', 'assoc_gene', 'tuning.rv', '.tsv', '.py', '/LARGE', 'nan',
            'NaN', '…', 'somatic', 'clonal'):
    if tok in ('nan', 'NaN') and re.search(rf'(?<![A-Za-z]){tok}(?![A-Za-z])', ms_text):
        print(f'    ✗ token in Multi-site carriers: {tok}'); bad += 1
for tok in ('_mainland', 'moderate_high', 'set_name', 'assoc_gene', 'tuning.rv', 'Control_', '.tsv', '.py',
            '/LARGE', 'nan', 'NaN', '…', 'somatic', 'clonal'):
    if re.search(rf'(?<![A-Za-z]){re.escape(tok)}(?![A-Za-z])' if tok in ('nan', 'NaN') else re.escape(tok), text):
        print(f'    ✗ token in the workbook: {tok}'); bad += 1
if ids & set(re.findall(r'[A-Za-z0-9_-]+', text)):
    print('    ✗ a sample ID outside sheet Multi-site carriers'); bad += 1
ms = pd.read_csv(f'{out}/gene_signals.multisite.tsv', sep='\t', dtype=str)
if set(ms.IID) - set(re.findall(r'[A-Za-z0-9_-]+', ms_text)):
    print('    ✗ a multi-site carrier missing from its sheet'); bad += 1
if (text + ms_text).count('Abbreviations:') < len(want) - 1:       # every sheet but Rare-variant QC
    print('    ✗ abbreviations missing on a sheet'); bad += 1
w = json.load(open(f'{out}/workbook.json'))
print(f"    {len(wb.sheetnames)} sheets; {w['variants']['rows']} variant rows; {len(ids)} sample IDs checked absent")
sys.exit(1 if bad else 0)
PY
ck "sheets, Excel-strict XML, no internal tokens / interpretation, sample IDs only in Multi-site carriers" $?

if [ $fail -eq 0 ]; then echo "═══ ALL PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
exit $fail
