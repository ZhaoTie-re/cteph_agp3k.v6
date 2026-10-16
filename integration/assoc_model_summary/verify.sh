#!/bin/bash
# Pre-flight and post-run checks for assoc_model_summary. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree and the upstream runs and must pass before
# a run. Sections 6-8 read results/ (or $OUT) and are skipped until it exists.
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/assoc_model_summary
cd "$C"
OUT=${OUT:-results}
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
source activate dsl2 2>/dev/null
out=$(nextflow lint assoc_model_summary.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
# the pipeline's own interpreter, first on PATH for every python3 below
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" assoc_model_summary.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import pandas, numpy, scipy, matplotlib, statsmodels" 2>/dev/null
ck "params.python ($PYBIN) imports pandas, numpy, scipy, matplotlib, statsmodels" $?

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' assoc_model_summary.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' assoc_model_summary.nf; ck "no task.cpus" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E 'clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M"; bad=1; }
done
ck "every clusterOptions m= is a multiple of 4571" $bad
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*(python3|Rscript|source activate)\s' assoc_model_summary.nf; ck "no bare python3 / Rscript / source activate in a process" $([ $? -eq 1 ] && echo 0 || echo 1)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('assoc_model_summary.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
cfg = pathlib.Path('nextflow.config').read_text()
bad = 0
for p, b in bodies.items():
    if p in ('INDEX', 'WRITE_RUN_MANIFEST'):
        ok = "executor 'local'" in b
    else:
        ok = "executor 'slurm'" in b and "queue 'gr10478b'" in b and re.search(r"\n\s*time '", b) \
            and re.search(rf"withName: '[^']*\b{p}\b[^']*'", cfg)
    if not ok:
        print(f'    ✗ {p}: executor / queue / time / resources'); bad += 1
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "executor/queue/time/resources per process" $?

echo "═══ 3. every script call's flags are accepted ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('assoc_model_summary.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
# input variable -> script file, from the workflow calls (script_file('x.py') in argument order)
SCRIPT_OF = {'script': None, 'plot_paired': 'plot_paired.py', 'plot_calibration': 'plot_calibration.py',
             'scan_tables': 'scan_tables.py', 'write_readme': 'write_readme.py', 'pc_effect_tables': 'pc_effect_tables.py',
             'plot_main_summary': 'plot_main_summary.py'}
SINGLE = {'REDUCE_SCAN': 'reduce_scan.py', 'PC_ASSOC': 'pc_assoc.py'}
bad, n = 0, 0
for proc, body in bodies.items():
    m = re.search(r'script:\s*"""(.*?)"""', body, re.S)
    if not m:
        continue
    calls = re.split(r'\$\{params\.python\}\s+\$\{(\w+)\}', m.group(1))
    for var, args in zip(calls[1::2], calls[2::2]):
        sc = SINGLE.get(proc) if var == 'script' else SCRIPT_OF.get(var)
        if not sc:
            print(f'  ✗ {proc}: unknown script variable {var}'); bad += 1; continue
        passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', args))
        h = subprocess.run([sys.executable, f'scripts/{sc}', '--help'], capture_output=True, text=True)
        known = set(re.findall(r'--([a-z][a-z0-9-]*)', h.stdout))
        miss = passed - known
        n += 1
        if h.returncode or miss:
            print(f'  ✗ {proc}: {sc} does not accept --' + ', --'.join(sorted(miss))); bad += 1
print(f'    {n} script calls checked')
sys.exit(1 if bad else 0)
PY
ck "flags accepted" $?

echo "═══ 4. upstream ═══"
P=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis
bad=0
for c in narrow_mainland intermediate_mainland full_mainland; do
  for f in $P/assoc_saige/results/$c/01.assoc/main/sumstats.tsv $P/assoc_saige/results/$c/01.assoc/detect/sumstats.tsv \
           $P/assoc_saige/results/$c/00.prep/pheno_covar.tsv $P/assoc_saige/results/$c/03.peaks/model_peaks_annotation.tsv \
           $P/assoc_plink2/results/$c/03.peaks/model_peaks_annotation.tsv $P/assoc_plink2/results/$c/02.scan/scan_qc.tsv; do
    [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }
  done
  n=$(find $P/assoc_plink2/results/$c/01.assoc/additive/ -maxdepth 1 -name '*.glm.logistic' -xtype f | wc -l)
  [ "$n" -eq 1 ] || { echo "    ✗ plink2 $c additive sumstats: $n files (symlinks into assoc_plink2/work?)"; bad=1; }
done
for f in $P/assoc_saige/results/_comparison/tables/lead_bymodel.tsv $P/assoc_saige/results/_comparison/tables/lead_crosscohort.tsv \
         $P/assoc_saige/results/_run_info/run_manifest.json $P/assoc_plink2/results/_run_info/run_manifest.json \
         $P/_shared/scripts/plot_style.py; do
  [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }
done
ck "upstream scans, peaks, tables, manifests and plot_style present" $bad

echo "═══ 5. docs ═══"
bad=0; for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md; do [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
ck "README, METHODS, OUTPUTS, STUDY_NOTES present" $bad

if [ ! -f "$OUT/_run_info/trace.txt" ]; then
  echo "═══ post-run sections skipped: no $OUT/_run_info/trace.txt ═══"
  if [ $fail -eq 0 ]; then echo "═══ PRE-FLIGHT PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
  exit $fail
fi

echo "═══ 6. DAG cardinality (trace.txt) ═══"
OUT="$OUT" python3 - <<'PY'
import csv, os, sys, collections
out = os.environ['OUT']
EXPECT = {'REDUCE_SCAN': 9, 'PC_ASSOC': 1, 'FIG_FIXED_RANDOM': 1, 'FIG_PC_EFFECT': 1, 'FIG_MAIN': 1, 'INDEX': 1,
          'WRITE_RUN_MANIFEST': 1}
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

echo "═══ 7. scan statistics against the upstream scan_qc ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, sys, glob
import pandas as pd
out = os.environ['OUT']
P = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis'
up = {'fixed': ('assoc_plink2', 'additive'), 'grm_pc': ('assoc_saige', 'main'), 'grm_only': ('assoc_saige', 'detect')}
bad, n = 0, 0
for f in sorted(glob.glob(f'{out}/00.scans/*.stats.json')):
    s = json.load(open(f))
    comp, model = up[s['model']]
    q = pd.read_csv(f'{P}/{comp}/results/_comparison/tables/scan_qc_all.tsv', sep='\t')
    q = q[(q.cohort == s['cohort']) & (q.model == model)].iloc[0]
    n += 1
    if abs(q.lambda_gc - s['lambda_gc']) > 1e-4 or int(q.n_suggestive) != s['n_p_lt_suggestive']:
        print(f"    ✗ {s['cohort']} {s['model']}: lambda {s['lambda_gc']:.5f} vs {q.lambda_gc}"); bad += 1
print(f'    {n} scans compared')
sys.exit(1 if bad or n != 9 else 0)
PY
ck "lambda_GC and P<1e-5 counts equal the upstream scan_qc" $?

echo "═══ 8. deliverables ═══"
bad=0
for f in 01.fixed_vs_random/fixed_vs_random.png 01.fixed_vs_random/calibration.png 02.pc_effect/grm_pc_vs_grm.png \
         03.main_results/main_manhattan.png 03.main_results/main_summary.png; do
  s=$(stat -c %s "$OUT/$f" 2>/dev/null || echo 0); [ "$s" -gt 100000 ] || { echo "    ✗ $f ($s bytes)"; bad=1; }
done
for f in README.md 01.fixed_vs_random/README.md 02.pc_effect/README.md 03.main_results/README.md; do
  [ -s "$OUT/$f" ] || { echo "    ✗ $f missing"; bad=1; continue; }
  grep -qwE 'nan|None|NaN' "$OUT/$f" && { echo "    ✗ $f contains nan/None"; bad=1; }
  [ "$f" = README.md ] || grep -q '^## Slide summary' "$OUT/$f" || { echo "    ✗ $f has no slide summary"; bad=1; }
done
grep -q 'Conclusion: pending' "$OUT/02.pc_effect/README.md" || { echo "    ✗ 02.pc_effect is no longer marked pending"; bad=1; }
ck "5 PNGs, 4 READMEs with slide summaries, no nan" $bad

if [ $fail -eq 0 ]; then echo "═══ ALL PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
exit $fail
